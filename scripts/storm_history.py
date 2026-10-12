#!/usr/bin/env python3
"""June-July wind history, 1976-2025, for each marine forecast area the race passes through.

Warning archives do not go back 50 years (Environment Canada keeps none public; the US National Weather Service's, via the Iowa
Environmental Mesonet, starts in 2005), so the 50-year record is rebuilt from wind: ERA5 reanalysis (ECMWF, via the Open-Meteo
archive API, hourly 10 m wind) at a point on the route in each area, counted against the warning thresholds:
  small craft: sustained 20-33 kn in Canada (Environment Canada strong wind warning), 23-33 kn in Alaska (NWS small craft advisory)
  gale 34-47 kn, storm 48-63 kn.
A day counts at the highest level any hour of it reaches. ERA5 is a ~30 km grid and smooths the wind in narrow straits, so where
a Fisheries and Oceans Canada (MEDS) weather buoy sits in the area, its June-July observations are counted the same way, and the
ERA5 hourly wind at the buoy is compared with it (ratio of their 95th percentiles) to give a corrected count as well.
For the US zone the actual small craft advisories, gale and storm warnings issued since 2005 are counted too.

Writes data/storm_history/{era5,buoys,nws,summary}.json and data/storm_history/area_year.csv.
"""
from __future__ import annotations

import csv
import io
import json
import math
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from datetime import date, datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "storm_history"
Y0, Y1 = 1976, 2025
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
# forecast area, country, a point on the route in it (lat, lon), the MEDS buoy in or next to it (if any), route note
AREAS = [
    ("Juan de Fuca Strait - east entrance", "CA", (48.36, -123.38), None, "Victoria start"),
    ("Haro Strait", "CA", (48.60, -123.22), None, ""),
    ("Strait of Georgia - south of Nanaimo", "CA", (49.25, -123.65), "46146", "Halibut Bank buoy"),
    ("Strait of Georgia - north of Nanaimo", "CA", (49.80, -124.75), "46131", "Sentry Shoal buoy"),
    ("Johnstone Strait", "CA", (50.42, -125.90), None, ""),
    ("Queen Charlotte Strait", "CA", (50.75, -127.20), None, ""),
    ("Queen Charlotte Sound - eastern half", "CA", (51.25, -127.95), "46204", "West Sea Otter buoy; Cape Caution crossing"),
    ("Central Coast from McInnes Island to Pine Island", "CA", (51.75, -128.00), None, "Fitz Hugh Sound to Bella Bella"),
    ("Inside channels, Finlayson to Grenville (no forecast area)", "CA", (53.00, -129.05), None, "inside route only; not a forecast area"),
    ("Hecate Strait - southern half", "CA", (52.70, -129.80), "46185", "South Hecate buoy; outside route only"),
    ("Hecate Strait - northern half", "CA", (53.70, -130.60), "46183", "North Hecate buoy; outside route only"),
    ("Dixon Entrance east", "CA", (54.55, -130.90), "46145", "Central Dixon Entrance buoy"),
    ("Clarence Strait (PKZ036)", "US", (55.10, -131.20), None, "Revillagigedo Channel to Ketchikan finish"),
]
# Buoy positions (approximate, for the ERA5 comparison)
BUOYS = {"46146": (49.34, -123.73), "46131": (49.91, -124.99), "46204": (51.38, -128.77), "46185": (52.42, -129.79),
         "46183": (53.62, -131.10), "46145": (54.37, -132.42)}
NWS_ZONES = {"Clarence Strait (PKZ036)": "PKZ036"}
# Observed wind: Environment Canada lightstations and automatic stations in each area (names as in its station inventory, in order of
# preference; each name's stations within 25 km of the area are joined into one record), and NDBC stations in Alaska.
OBS = {
    "Juan de Fuca Strait - east entrance": ["VICTORIA GONZALES", "TRIAL ISLAND"],
    "Haro Strait": ["KELP REEFS", "DISCOVERY ISLAND"],
    "Strait of Georgia - south of Nanaimo": ["ENTRANCE ISLAND", "BALLENAS ISLAND"],
    "Strait of Georgia - north of Nanaimo": ["SISTERS ISLAND", "GRIEF POINT", "MERRY ISLAND"],
    "Johnstone Strait": ["FANNY ISLAND", "HELMCKEN ISLAND", "CHATHAM POINT"],
    "Queen Charlotte Strait": ["HERBERT ISLAND", "PINE ISLAND", "ALERT BAY"],
    "Queen Charlotte Sound - eastern half": ["EGG ISLAND"],
    "Central Coast from McInnes Island to Pine Island": ["ADDENBROKE ISLAND", "DRYAD POINT", "BELLA BELLA", "IVORY ISLAND", "MCINNES ISLAND"],
    "Inside channels, Finlayson to Grenville (no forecast area)": ["BOAT BLUFF"],
    "Hecate Strait - southern half": [],
    "Hecate Strait - northern half": ["BONILLA ISLAND"],
    "Dixon Entrance east": ["GREY ISLET", "GREEN ISLAND", "TRIPLE ISLAND", "LUCY ISLAND", "HOLLAND ROCK"],
    "Clarence Strait (PKZ036)": ["ndbc:mrya2", "ndbc:ktna2"],
}
OBS_KM = 80   # how far from the area's route point a named station may be


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, timeout=45, tries=8):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read()
        except Exception as e:
            log(f"  {type(e).__name__} {e} ({url[:90]})")
            time.sleep(min(120, 10 * (i + 1) ** 2))
    return None


def thresholds(country):
    return (23 if country == "US" else 20), 34, 48


def level(kn, country):
    sc, g, s = thresholds(country)
    return 3 if kn >= s else 2 if kn >= g else 1 if kn >= sc else 0


HOURS = 61 * 24   # June 1 00:00 to July 31 24:00, Pacific daylight time (UTC-7)
MERGE_H, MIN_H = 6, 2   # hours at small craft or more, with gaps up to MERGE_H, make one event; events shorter than MIN_H are dropped


def era5(lat, lon, year):
    """Hourly 10 m wind, kn, June 1 - July 31 of `year`, a list of HOURS values (None where missing)."""
    q = urllib.parse.urlencode({"latitude": lat, "longitude": lon, "start_date": f"{year}-06-01", "end_date": f"{year}-07-31",
                                "hourly": "wind_speed_10m", "wind_speed_unit": "kn", "models": "era5", "timezone": "Etc/GMT+7"})
    b = get("https://archive-api.open-meteo.com/v1/archive?" + q)
    if not b:
        return None
    v = json.loads(b)["hourly"]["wind_speed_10m"]
    return [None if x is None else round(x, 1) for x in v[:HOURS]] + [None] * max(0, HOURS - len(v))


def events(hourly, country, scale=1.0):
    """[[start hour, hours, level]] of runs at small craft or more (level = highest reached)."""
    out, cur, gap = [], None, 0
    for h, x in enumerate(hourly):
        lv = level(x * scale, country) if x is not None else 0
        if lv:
            if cur and h - cur[1] <= MERGE_H:
                cur[1] = h; cur[2] = max(cur[2], lv)
            else:
                if cur: out.append(cur)
                cur = [h, h, lv]
    if cur: out.append(cur)
    return [[a, b - a + 1, lv] for a, b, lv in out if b - a + 1 >= MIN_H]


def pct(vals, p):
    v = sorted(vals)
    return v[min(len(v) - 1, int(p / 100 * len(v)))] if v else None


def meds(buoy):
    """MEDS buoy hourly wind (kn) June-July, {year: list of HOURS}; WSPD is m/s (the higher of the two anemometers is used)."""
    for url in (f"https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{buoy}_csv.zip",
                f"https://meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{buoy}_csv.zip"):
        b = get(url, 300, 3)
        if b:
            break
    else:
        return None
    z = zipfile.ZipFile(io.BytesIO(b))
    out = {}
    with z.open(z.namelist()[0]) as f:
        rd = csv.DictReader(io.TextIOWrapper(f, "latin-1"))
        for r in rd:
            t = (r.get("DATE") or r.get("Date") or "").strip()
            try:
                dt = datetime.strptime(t[:16], "%m/%d/%Y %H:%M")
            except ValueError:
                try:
                    dt = datetime.fromisoformat(t[:16])
                except ValueError:
                    continue
            dt -= timedelta(hours=7)   # UTC to PDT, as ERA5 above
            if dt.month not in (6, 7) or r.get("Q_FLAG", "1") in ("3", "4", "9"):
                continue
            w = []
            for k in ("WSPD", "WSS$", "WSPD$", "WSS"):
                try:
                    x = float(r.get(k, ""))
                    if 0 <= x < 60:
                        w.append(x)
                except ValueError:
                    pass
            if not w:
                continue
            h = int((dt - datetime(dt.year, 6, 1)).total_seconds() // 3600)
            if 0 <= h < HOURS:
                out.setdefault(dt.year, [None] * HOURS)[h] = round(max(w) * 1.943844, 1)
    log(f"  buoy {buoy}: {len(out)} years, header {rd.fieldnames[:14]}")
    return out


def nws(zone):
    """NWS small craft advisories and gale and storm warnings issued for a zone since 2005: June-July events as [[year, start hour, hours, level]]."""
    b = get(f"https://mesonet.agron.iastate.edu/json/vtec_events_byugc.php?ugc={zone}&sdate=2005-01-01&edate={Y1}-12-31")
    if not b:
        return None
    out = []
    for e in json.loads(b).get("events", []):
        lv = {("SC", "Y"): 1, ("GL", "W"): 2, ("SR", "W"): 3, ("HF", "W"): 3}.get((e.get("phenomena"), e.get("significance")))
        if not lv:
            continue
        i = datetime.fromisoformat(e["issue"][:19]) - timedelta(hours=7)
        x = datetime.fromisoformat(e["expire"][:19]) - timedelta(hours=7)
        h0 = int((i - datetime(i.year, 6, 1)).total_seconds() // 3600); h1 = int((x - datetime(i.year, 6, 1)).total_seconds() // 3600)
        if h1 < 0 or h0 >= HOURS:
            continue
        out.append([i.year, max(0, h0), min(HOURS, h1) - max(0, h0), lv])
    return out


INV = "https://collaboration.cmc.ec.gc.ca/cmc/climate/Get_More_Data_Plus_de_donnees/Station%20Inventory%20EN.csv"


def inventory():
    b = get(INV, 120, 4)
    if not b:
        return []
    lines = b.decode("utf-8-sig", "replace").splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith('"Name"'))
    return list(csv.DictReader(lines[k:]))


def km(a, b):
    return nm(a, b) * 1.852


def nm(a, b):
    return math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0])))


def fill(series):
    """Carry each observation forward up to 3 h, so 3- and 4-hourly lightstation reports read as a continuous record."""
    out, last, age = list(series), None, 9
    for h, x in enumerate(series):
        if x is not None:
            last, age = x, 0
        else:
            age += 1
            if last is not None and age <= 3:
                out[h] = last
    return out


def eccc(st, y0, y1):
    """{year: list of HOURS kn (PDT)} from one ECCC station's hourly data, June and July."""
    out = {}
    sid = st["Station ID"]
    for y in range(max(Y0, y0), min(Y1, y1) + 1):
        ser = [None] * HOURS
        for m in (6, 7):
            b = get(f"https://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID={sid}&Year={y}&Month={m}&Day=1"
                    f"&time=LST&timeframe=1&submit=Download+Data", 90, 4)
            if not b:
                continue
            for r in csv.DictReader(io.StringIO(b.decode("utf-8-sig", "replace"))):
                v = r.get("Wind Spd (km/h)") or ""
                try:
                    kn = float(v) / 1.852
                    dt = datetime.strptime(r["Date/Time (LST)"][:16], "%Y-%m-%d %H:%M") + timedelta(hours=1)   # LST (UTC-8) to PDT
                except (ValueError, KeyError, TypeError):
                    continue
                h = int((dt - datetime(y, 6, 1)).total_seconds() // 3600)
                if 0 <= h < HOURS:
                    ser[h] = round(kn, 1)
            time.sleep(0.3)
        if sum(x is not None for x in ser) > 24 * 10:
            out[y] = ser
    return out


def ndbc(sid):
    """{year: list of HOURS kn (PDT)} from an NDBC station's historical standard meteorological files."""
    import gzip
    out = {}
    for y in range(Y0, Y1 + 1):
        b = get(f"https://www.ndbc.noaa.gov/data/historical/stdmet/{sid}h{y}.txt.gz", 60, 2)
        if not b:
            continue
        try:
            txt = gzip.decompress(b).decode("latin-1")
        except OSError:
            continue
        ser = [None] * HOURS
        for line in txt.splitlines():
            p = line.split()
            if not p or p[0].startswith("#") or not p[0].isdigit():
                continue
            try:
                yy = int(p[0]); yy += 1900 if yy < 100 else 0
                dt = datetime(yy, int(p[1]), int(p[2]), int(p[3])) - timedelta(hours=7)
                w = float(p[6] if len(p) > 16 else p[5])   # WSPD m/s (column 6 once minutes were added in 2005)
            except (ValueError, IndexError):
                continue
            if w >= 90:
                continue
            h = int((dt - datetime(y, 6, 1)).total_seconds() // 3600)
            if 0 <= h < HOURS:
                ser[h] = round(max(ser[h] or 0, w * 1.943844), 1)
        if sum(x is not None for x in ser) > 24 * 10:
            out[y] = ser
        time.sleep(0.3)
    return out


def observed(name, point):
    """The area's observed June-July record: {year: hourly kn}, joined from its stations in order of preference, with the source per year."""
    inv = observed.inv if hasattr(observed, "inv") else None
    rec, src = {}, {}
    for want in OBS.get(name, []):
        if want.startswith("ndbc:"):
            parts = {want[5:].upper(): ndbc(want[5:])}
        else:
            if inv is None:
                inv = observed.inv = inventory()
            parts = {}
            for st in inv:
                try:
                    ll = (float(st["Latitude (Decimal Degrees)"]), float(st["Longitude (Decimal Degrees)"]))
                    h0, h1 = int(st["HLY First Year"] or 0), int(st["HLY Last Year"] or 0)
                except (ValueError, KeyError):
                    continue
                if want in st["Name"].upper() and st["Province"].startswith("BRITISH") and h1 >= Y0 and km(ll, point) < OBS_KM:
                    log(f"  {name}: {st['Name']} ({st['Station ID']}) hourly {h0}-{h1}, {km(ll, point):.0f} km")
                    parts[st["Name"]] = eccc(st, h0, h1)
        for sname, ser in parts.items():
            for y, v in ser.items():
                v = fill(v)
                if sum(x is not None for x in v) >= 0.4 * HOURS and y not in rec:
                    rec[y] = v; src[y] = sname
    return rec, src


def stats(ev_by_year, hourly_by_year, country, scale):
    """Likelihood, frequency and duration of events at small craft or more (sc) and gale or more (gale)."""
    yrs = sorted(ev_by_year)
    n = len(yrs)
    s = {}
    for key, lv in (("sc", 1), ("gale", 2), ("storm", 3)):
        evs = [(y, e) for y in yrs for e in ev_by_year[y] if e[2] >= lv]
        dur = [e[1] for _, e in evs]
        # chance of at least one such event starting in a 16-day stretch, over every 16-day stretch of June-July
        hit = tot = 0
        for y in yrs:
            st = [e[0] for e in ev_by_year[y] if e[2] >= lv]
            for d0 in range(0, 61 - 16 + 1):
                tot += 1; hit += any(d0 * 24 <= h < (d0 + 16) * 24 for h in st)
        days = [len({(e[0] + k) // 24 for e in ev_by_year[y] if e[2] >= lv for k in range(e[1])}) for y in yrs]
        hrs = sum(sum(1 for x in hourly_by_year[y] if x is not None and level(x * scale, country) >= lv) for y in yrs)
        s[key] = {"events_per_jj": round(len(evs) / n, 2), "days_per_jj": round(sum(days) / n, 1), "share_of_hours_pct": round(100 * hrs / (n * HOURS), 2),
                  "p_in_16_days": round(hit / tot, 3) if tot else None, "years_with_any_pct": round(100 * sum(d > 0 for d in days) / n),
                  "dur_h": {"p25": pct(dur, 25), "median": pct(dur, 50), "p75": pct(dur, 75), "p90": pct(dur, 90), "max": max(dur) if dur else None},
                  "longest_by_days": {"<12 h": sum(d < 12 for d in dur), "12-24 h": sum(12 <= d < 24 for d in dur), "1-2 days": sum(24 <= d < 48 for d in dur),
                                      "2-3 days": sum(48 <= d < 72 for d in dur), "3+ days": sum(d >= 72 for d in dur)},
                  "by_decade_events": {f"{d}-{d + 9}": round(sum(1 for y, _ in evs if d <= y < d + 10) / max(1, sum(1 for y in yrs if d <= y < d + 10)), 2) for d in range(Y0, Y1 + 1, 10)},
                  "by_half_month_share_pct": {lab: round(100 * sum(1 for y in yrs for h in range(a * 24, b * 24) if hourly_by_year[y][h] is not None and level(hourly_by_year[y][h] * scale, country) >= lv)
                                                        / (n * (b - a) * 24), 2) for lab, a, b in (("Jun 1-15", 0, 15), ("Jun 16-30", 15, 30), ("Jul 1-15", 30, 45), ("Jul 16-31", 45, 61))}}
    return s


def points():
    return [p for name, cc, (la, lo), buoy, note in AREAS for p in [(name, la, lo)] + ([(f"buoy {buoy}", *BUOYS[buoy])] if buoy else [])]


def main():
    """--only N: fetch ERA5 for the Nth point alone into parts/N.json (one runner each: Open-Meteo throttles a single address);
    otherwise merge any parts into the cache, fetch whatever is still missing, and work out the history."""
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "era5_hourly.json"
    raw = json.loads(cache.read_text()) if cache.exists() else {}
    parts = OUT / "parts"
    if "--list" in sys.argv:   # stations with hourly data near each area, to choose OBS from
        inv = inventory()
        for name, cc, pt, _, _ in AREAS:
            log(f"== {name}")
            for st in sorted(inv, key=lambda st: km(pt, (float(st["Latitude (Decimal Degrees)"] or 0), float(st["Longitude (Decimal Degrees)"] or 0)))):
                ll = (float(st["Latitude (Decimal Degrees)"] or 0), float(st["Longitude (Decimal Degrees)"] or 0))
                if km(pt, ll) > 80:
                    break
                if st.get("HLY First Year"):
                    log(f"   {st['Name']:40} {st['Station ID']:>6} {ll[0]:.3f},{ll[1]:.3f} {km(pt, ll):5.0f} km  HLY {st['HLY First Year']}-{st['HLY Last Year']}")
        return
    if "--obs" in sys.argv:
        n = int(sys.argv[sys.argv.index("--obs") + 1]); name, cc, pt, _, _ = AREAS[n]
        rec, src = observed(name, pt)
        parts.mkdir(exist_ok=True)
        (parts / f"obs{n}.json").write_text(json.dumps({name: {"years": {str(y): v for y, v in rec.items()}, "src": {str(y): s for y, s in src.items()}}}, separators=(",", ":")))
        log(f"observed {name}: {len(rec)} years {sorted(set(src.values()))}")
        return
    if "--only" in sys.argv:
        n = int(sys.argv[sys.argv.index("--only") + 1]); key, a, b = points()[n]
        got = raw.get(key, {})
        for y in range(Y0, Y1 + 1):
            if str(y) not in got:
                d = era5(a, b, y)
                if d is not None:
                    got[str(y)] = d
                time.sleep(1)
        parts.mkdir(exist_ok=True)
        (parts / f"{n}.json").write_text(json.dumps({key: got}, separators=(",", ":")))
        log(f"ERA5 {key}: {len(got)} years")
        return
    ocache = OUT / "obs_hourly.json"
    obs = json.loads(ocache.read_text()) if ocache.exists() else {}
    for f in sorted(parts.glob("[0-9]*.json")) if parts.exists() else []:
        for k, v in json.loads(f.read_text()).items():
            raw.setdefault(k, {}).update(v)
    for f in sorted(parts.glob("obs*.json")) if parts.exists() else []:
        for k, v in json.loads(f.read_text()).items():
            if v["years"]:
                obs[k] = v
    cache.write_text(json.dumps(raw, separators=(",", ":")))
    ocache.write_text(json.dumps(obs, separators=(",", ":")))
    for key, a, b in points():   # anything a part did not bring
        raw.setdefault(key, {})
        for y in range(Y0, Y1 + 1):
            if str(y) in raw[key]:
                continue
            d = era5(a, b, y)
            if d is None:
                log(f"  no ERA5 for {key} {y}"); continue
            raw[key][str(y)] = d
            time.sleep(2.5)
        cache.write_text(json.dumps(raw, separators=(",", ":")))
        log(f"ERA5 {key}: {len(raw[key])} years")

    summary, out_ev, rows = [], [], []
    for ai, (name, cc, (la, lo), buoy, note) in enumerate(AREAS):
        hourly = {int(y): v for y, v in raw.get(name, {}).items()}
        if not hourly:
            continue
        scale, cal = 1.0, None
        bd = meds(buoy) if buoy else None
        if bd:
            era_b = {int(y): v for y, v in raw.get(f"buoy {buoy}", {}).items()}
            ob, md = [], []
            for y, bser in bd.items():
                if y not in era_b:
                    continue
                for d in range(61):
                    o = [x for x in bser[d * 24:(d + 1) * 24] if x is not None]; m = [x for x in era_b[y][d * 24:(d + 1) * 24] if x is not None]
                    if len(o) >= 18 and len(m) == 24:
                        ob.append(max(o)); md.append(max(m))
            if len(ob) > 300:
                scale = pct(ob, 95) / pct(md, 95)
                cal = {"buoy": buoy, "years": f"{min(bd)}-{max(bd)}", "days_compared": len(ob), "scale": round(scale, 3),
                       "buoy_p95_daily_max_kn": pct(ob, 95), "era5_p95_daily_max_kn": pct(md, 95),
                       "buoy_sc_days_pct": round(100 * sum(level(x, cc) >= 1 for x in ob) / len(ob), 1),
                       "era5_sc_days_pct": round(100 * sum(level(x, cc) >= 1 for x in md) / len(md), 1),
                       "era5_scaled_sc_days_pct": round(100 * sum(level(x * scale, cc) >= 1 for x in md) / len(md), 1)}
        if not 0.8 <= scale <= 1.25:   # ERA5 far off the buoy here (its cell is part land): not worth scaling
            scale = 1.0
        ev = {y: events(hourly[y], cc, scale) for y in sorted(hourly)}
        raw_ev = {y: events(hourly[y], cc) for y in sorted(hourly)}
        s = {"area": name, "country": cc, "point": [la, lo], "note": note, "years": len(ev), "scale": round(scale, 3), "calibration": cal,
             "thresholds_kn": dict(zip(("small_craft", "gale", "storm"), thresholds(cc))),
             "corrected": stats(ev, hourly, cc, scale), "era5_raw": stats(raw_ev, hourly, cc, 1.0)}
        if bd:
            bev = {y: events(v, cc) for y, v in bd.items() if sum(x is not None for x in v) > 0.75 * HOURS}
            if bev:
                s["buoy_observed"] = stats(bev, {y: bd[y] for y in bev}, cc, 1.0)
        if name in NWS_ZONES:
            nw = nws(NWS_ZONES[name])
            if nw is not None:
                byy = {y: [[h, d, lv] for yy, h, d, lv in nw if yy == y] for y in range(2005, Y1 + 1)}
                s["nws_issued_2005_on"] = stats(byy, {y: [None] * HOURS for y in byy}, cc, 1.0)
        # the record used for the area: the NWS warnings issued (Alaska), else its stations' observations (unless the buoy has the longer
        # record and the stations under 15 years), else the buoy, else ERA5
        ob = obs.get(name)
        if ob and len(ob["years"]) >= 8:
            oh = {int(y): v for y, v in ob["years"].items()}
            oev = {y: events(v, cc) for y, v in oh.items()}
            st = sorted(set(ob["src"].values()))
            s["observed"] = {**stats(oev, oh, cc, 1.0), "stations": st, "years": len(oh), "first": min(oh), "last": max(oh)}
        bev = bev if bd and "buoy_observed" in s else {}
        if "nws_issued_2005_on" in s:
            # in Alaska the warnings actually issued for the zone: the one NDBC station near the route reads far lighter
            best, bsrc, bh = byy, "NWS warnings issued", {y: [None] * HOURS for y in byy}
        elif "observed" in s and (len(oh) >= 15 or len(bev) < len(oh)):
            best, bsrc, bh = oev, "observed", oh
        elif len(bev) >= 8:
            best, bsrc, bh = bev, f"buoy {buoy}", {y: bd[y] for y in bev}
        else:
            best, bsrc, bh = ev, "ERA5" + (" (scaled)" if scale != 1.0 else ""), hourly
        s["best_source"] = bsrc; s["best_years"] = f"{min(best)}-{max(best)} ({len(best)} years)"
        s["best"] = s["observed"] if bsrc == "observed" else s.get("buoy_observed") if bsrc.startswith("buoy") \
            else s["nws_issued_2005_on"] if bsrc.startswith("NWS") else s["corrected"]
        summary.append(s)
        out_ev.append({"area": name, "country": cc, "point": [la, lo], "note": note, "source": bsrc, "years": {str(y): e for y, e in best.items()}})
        for y in sorted(set(best) | set(ev)):
            r = {"area": name, "year": y}
            for tag, E, Hs in (("best", best, bh), ("era5", ev, hourly)):
                if y in E:
                    r.update({f"{tag}_{k}_events": sum(e[2] >= lv for e in E[y]) for k, lv in (("sc", 1), ("gale", 2), ("storm", 3))})
                    r.update({f"{tag}_sc_hours": sum(e[1] for e in E[y]), f"{tag}_longest_h": max((e[1] for e in E[y]), default=0),
                              f"{tag}_max_kn": round(max((x for x in Hs[y] if x is not None), default=0) * (scale if tag == "era5" else 1), 1)})
            r["best_source"] = (ob["src"].get(str(y), "") if bsrc == "observed" else bsrc) if y in best else ""
            rows.append(r)
        log(name, bsrc, json.dumps({k: (s["best"][k]["events_per_jj"], s["best"][k]["dur_h"]["median"]) for k in ("sc", "gale", "storm")}))
    with open(OUT / "area_year.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(dict.fromkeys(k for r in rows for k in r))); w.writeheader(); w.writerows(rows)
    gen = time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime())
    (OUT / "summary.json").write_text(json.dumps({
        "note": "June 1 - July 31, 1976-2025, per marine forecast area on the route. best: the record used, best_source says which: observed (Environment "
                "Canada lightstation and automatic station hourly wind, NDBC in Alaska; years at least 40% complete once each report is carried 3 h; at least 8 years), else a MEDS buoy, else the NWS warnings issued (Clarence Strait, 2005 on), else ERA5 hourly 10 m wind "
                "at a route point in the area (era5_raw; corrected = scaled to the area's buoy where that is within 20%). An event is a run "
                f"of hours at small craft strength or more (gaps up to {MERGE_H} h merged, at least {MIN_H} h); its level is the highest it reached. "
                "corrected: ERA5 scaled by the buoy/ERA5 ratio of 95th percentile daily maximum wind where a MEDS buoy is in the area (else as era5_raw). "
                "buoy_observed: the buoy itself. nws_issued_2005_on: the warnings the US National Weather Service actually issued. "
                "Thresholds: small craft 20 kn (Canada) / 23 kn (US), gale 34 kn, storm 48 kn. Times Pacific daylight time.",
        "generated": gen, "areas": summary}, indent=1))
    (OUT / "events.json").write_text(json.dumps({
        "note": "Per area and year, June-July events at small craft strength or more: [start hour after June 1 00:00 PDT, hours, level 1 small craft / 2 gale / 3 storm]. "
                "From scripts/storm_history.py: each area's best record (source).", "generated": gen, "y0": Y0, "y1": Y1, "areas": out_ev}, separators=(",", ":")))
    log("done")


if __name__ == "__main__":
    main()
