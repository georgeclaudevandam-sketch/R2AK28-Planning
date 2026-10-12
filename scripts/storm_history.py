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
    ("Hecate Strait - southern half", "CA", (52.70, -129.80), "46185", "South Hecate buoy; outside route only"),
    ("Hecate Strait - northern half", "CA", (53.70, -130.60), "46183", "North Hecate buoy; outside route; nearest area to Grenville Channel"),
    ("Dixon Entrance east", "CA", (54.55, -130.90), "46145", "Central Dixon Entrance buoy"),
    ("Clarence Strait (PKZ036)", "US", (55.10, -131.20), None, "Revillagigedo Channel to Ketchikan finish"),
]
# Buoy positions (approximate, for the ERA5 comparison)
BUOYS = {"46146": (49.34, -123.73), "46131": (49.91, -124.99), "46204": (51.38, -128.77), "46185": (52.42, -129.79),
         "46183": (53.62, -131.10), "46145": (54.37, -132.42)}
NWS_ZONES = {"Clarence Strait (PKZ036)": "PKZ036"}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, timeout=120, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read()
        except Exception as e:
            log(f"  {type(e).__name__} {e} ({url[:90]})")
            time.sleep(min(300, 20 * (i + 1) ** 2))
    return None


def thresholds(country):
    return (23 if country == "US" else 20), 34, 48


def level(kn, country):
    sc, g, s = thresholds(country)
    return 3 if kn >= s else 2 if kn >= g else 1 if kn >= sc else 0


def era5(lat, lon, year):
    """Hourly 10 m wind, kn, June 1 - July 31 of `year`, as {date: [hourly kn]} (local standard time days)."""
    q = urllib.parse.urlencode({"latitude": lat, "longitude": lon, "start_date": f"{year}-06-01", "end_date": f"{year}-07-31",
                                "hourly": "wind_speed_10m", "wind_speed_unit": "kn", "models": "era5", "timezone": "America/Vancouver"})
    b = get("https://archive-api.open-meteo.com/v1/archive?" + q)
    if not b:
        return None
    h = json.loads(b)["hourly"]
    days = {}
    for t, v in zip(h["time"], h["wind_speed_10m"]):
        if v is not None:
            days.setdefault(t[:10], []).append(v)
    return days


def count(days, country, scale=1.0):
    """Days at small craft or higher, gale or higher, storm or higher; and episodes (runs of consecutive days at small craft+)."""
    lv = {d: level(max(v) * scale, country) for d, v in days.items() if v}
    n = [sum(1 for x in lv.values() if x >= k) for k in (1, 2, 3)]
    ep, prev = 0, False
    d0 = min(lv) if lv else None
    for d in sorted(lv):
        on = lv[d] >= 1
        ep += on and not prev
        prev = on
    return {"sc_days": n[0], "gale_days": n[1], "storm_days": n[2], "episodes": ep, "days": len(lv), "max_kn": round(max((max(v) for v in days.values()), default=0) * scale, 1)}


def pct(vals, p):
    v = sorted(vals)
    return v[min(len(v) - 1, int(p / 100 * len(v)))] if v else None


def meds(buoy):
    """MEDS buoy hourly wind (kn) June-July, {date: [kn]}; WSPD is m/s (the higher of the two anemometers is used)."""
    for url in (f"https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{buoy}_csv.zip",
                f"https://meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{buoy}_csv.zip"):
        b = get(url, 300, 3)
        if b:
            break
    else:
        return None
    z = zipfile.ZipFile(io.BytesIO(b))
    days = {}
    with z.open(z.namelist()[0]) as f:
        rd = csv.DictReader(io.TextIOWrapper(f, "latin-1"))
        for r in rd:
            t = r.get("DATE") or r.get("Date") or ""
            try:
                dt = datetime.strptime(t.strip()[:16], "%m/%d/%Y %H:%M")
            except ValueError:
                try:
                    dt = datetime.fromisoformat(t.strip()[:16])
                except ValueError:
                    continue
            dt -= timedelta(hours=8)   # UTC to PST days, as ERA5 above
            if dt.month not in (6, 7):
                continue
            w = []
            for k in ("WSPD", "WSS$", "WSPD$", "WSS"):
                try:
                    x = float(r.get(k, ""))
                    if 0 <= x < 60:
                        w.append(x)
                except ValueError:
                    pass
            if r.get("Q_FLAG", "0") not in ("0", "1", "", None):
                continue
            if w:
                days.setdefault(dt.date().isoformat(), []).append(max(w) * 1.943844)
    log(f"  buoy {buoy}: {len(days)} June-July days, header {rd.fieldnames[:12]}")
    return days


def nws(zone):
    b = get(f"https://mesonet.agron.iastate.edu/json/vtec_events_byugc.php?ugc={zone}&sdate=2005-01-01&edate={Y1}-12-31")
    if not b:
        return None
    ev = json.loads(b).get("events", [])
    out = {}
    for e in ev:
        ph, sg = e.get("phenomena"), e.get("significance")
        kind = {("SC", "Y"): "sc", ("GL", "W"): "gale", ("SR", "W"): "storm", ("HF", "W"): "storm"}.get((ph, sg))
        if not kind:
            continue
        i = datetime.fromisoformat(e["issue"].replace("Z", "+00:00")[:19]) - timedelta(hours=8)
        x = datetime.fromisoformat(e["expire"].replace("Z", "+00:00")[:19]) - timedelta(hours=8)
        d = i.date()
        while d <= x.date():
            if d.month in (6, 7):
                rank = {"sc": 1, "gale": 2, "storm": 3}[kind]
                out[d.isoformat()] = max(out.get(d.isoformat(), 0), rank)
            d += timedelta(days=1)
    return {"days": out, "events": len(ev)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "era5_daily_max.json"
    raw = json.loads(cache.read_text()) if cache.exists() else {}
    for name, cc, (la, lo), buoy, note in AREAS:
        pts = [(name, la, lo)] + ([(f"buoy {buoy}", *BUOYS[buoy])] if buoy else [])
        for key, a, b in pts:
            raw.setdefault(key, {})
            for y in range(Y0, Y1 + 1):
                if str(y) in raw[key]:
                    continue
                d = era5(a, b, y)
                if d is None:
                    log(f"  no ERA5 for {key} {y}"); continue
                raw[key][str(y)] = {k: [round(max(v), 1), round(sorted(v)[len(v) // 2], 1)] for k, v in d.items()}   # daily max, median
                time.sleep(2.5)
            cache.write_text(json.dumps(raw, separators=(",", ":")))
            log(f"ERA5 {key}: {len(raw[key])} years")

    def days_of(key, years=None):
        return {d: [v[0]] for y, dd in raw.get(key, {}).items() if years is None or int(y) in years for d, v in dd.items()}

    buoys = {}
    for name, cc, _, buoy, _ in AREAS:
        if buoy:
            bd = meds(buoy)
            if bd:
                buoys[buoy] = bd
    rows, summary = [], []
    for name, cc, (la, lo), buoy, note in AREAS:
        scale, cal = 1.0, None
        if buoy in buoys:
            bd = buoys[buoy]
            yrs = sorted({int(d[:4]) for d in bd})
            era_b = days_of(f"buoy {buoy}", set(yrs))
            both = [d for d in bd if d in era_b and len(bd[d]) >= 12]
            ob = [max(bd[d]) for d in both]; md = [era_b[d][0] for d in both]
            if len(both) > 200:
                scale = pct(ob, 95) / pct(md, 95)
                bo = count({d: bd[d] for d in both}, cc); em = count({d: era_b[d] for d in both}, cc)
                cal = {"buoy": buoy, "years": f"{yrs[0]}-{yrs[-1]}", "days_compared": len(both), "scale": round(scale, 2),
                       "buoy_sc_days_per_jj": round(bo["sc_days"] / len(both) * 61, 1), "era5_sc_days_per_jj": round(em["sc_days"] / len(both) * 61, 1),
                       "buoy_gale_days_per_jj": round(bo["gale_days"] / len(both) * 61, 2), "era5_gale_days_per_jj": round(em["gale_days"] / len(both) * 61, 2)}
        nw = nws(NWS_ZONES[name]) if name in NWS_ZONES else None
        per = []
        for y in range(Y0, Y1 + 1):
            dd = {d: [v[0]] for d, v in raw.get(name, {}).get(str(y), {}).items()}
            if not dd:
                continue
            r = count(dd, cc); c = count(dd, cc, scale)
            row = {"area": name, "year": y, **{f"era5_{k}": v for k, v in r.items() if k != "days"},
                   **{f"corr_{k}": v for k, v in c.items() if k in ("sc_days", "gale_days", "storm_days", "episodes")}}
            if buoy in buoys:
                bdy = {d: v for d, v in buoys[buoy].items() if d.startswith(str(y)) and len(v) >= 12}
                if len(bdy) >= 45:
                    b = count(bdy, cc)
                    row.update(buoy_days=len(bdy), buoy_sc_days=b["sc_days"], buoy_gale_days=b["gale_days"], buoy_storm_days=b["storm_days"], buoy_max_kn=b["max_kn"])
            if nw and y >= 2005:
                yd = [v for d, v in nw["days"].items() if d.startswith(str(y))]
                row.update(nws_sc_days=sum(v >= 1 for v in yd), nws_gale_days=sum(v >= 2 for v in yd), nws_storm_days=sum(v >= 3 for v in yd))
            rows.append(row); per.append(row)
        if not per:
            continue
        n = len(per)
        avg = lambda k: round(sum(r.get(k, 0) for r in per) / n, 2)
        dec = {}
        for d0 in range(1976, 2026, 10):
            sub = [r for r in per if d0 <= r["year"] < d0 + 10]
            if sub:
                dec[f"{d0}s"] = round(sum(r["corr_sc_days"] for r in sub) / len(sub), 1)
        # chance a given June-July day is at small craft or higher, by half-month
        halves = {}
        for lab, m, a, b in (("Jun 1-15", 6, 1, 15), ("Jun 16-30", 6, 16, 30), ("Jul 1-15", 7, 1, 15), ("Jul 16-31", 7, 16, 31)):
            tot = hit = 0
            for y, dd in raw.get(name, {}).items():
                for d, v in dd.items():
                    if int(d[5:7]) == m and a <= int(d[8:]) <= b:
                        tot += 1; hit += level(v[0] * scale, cc) >= 1
            halves[lab] = round(hit / tot, 3) if tot else None
        s = {"area": name, "country": cc, "point": [la, lo], "note": note, "years": n,
             "era5_sc_days": avg("era5_sc_days"), "era5_gale_days": avg("era5_gale_days"), "era5_storm_days": avg("era5_storm_days"),
             "corr_sc_days": avg("corr_sc_days"), "corr_gale_days": avg("corr_gale_days"), "corr_storm_days": avg("corr_storm_days"),
             "corr_episodes": avg("corr_episodes"), "years_with_gale_pct": round(100 * sum(r["corr_gale_days"] > 0 for r in per) / n),
             "max_kn_era5": max(r["era5_max_kn"] for r in per), "by_decade_corr_sc_days": dec, "day_chance_sc_corr": halves, "calibration": cal}
        if nw:
            sub = [r for r in per if r["year"] >= 2005]
            s["nws_2005_2025"] = {k: round(sum(r.get(k, 0) for r in sub) / len(sub), 2) for k in ("nws_sc_days", "nws_gale_days", "nws_storm_days")}
            s["nws_events_all_months"] = nw["events"]
        summary.append(s)
        log(name, json.dumps({k: s[k] for k in ("corr_sc_days", "corr_gale_days", "corr_storm_days")}), "cal", cal)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("area", "year"), k))
    keys = ["area", "year"] + [k for k in keys if k not in ("area", "year")]
    with open(OUT / "area_year.csv", "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    (OUT / "summary.json").write_text(json.dumps({
        "note": "June 1 - July 31, 1976-2025. Days per June-July at or above each warning level from ERA5 hourly 10 m wind at a route point "
                "per marine forecast area (era5_*), the same scaled by the buoy/ERA5 ratio of 95th percentile daily max wind where a MEDS "
                "buoy is in the area (corr_*; equal to era5_* elsewhere), MEDS buoy observations (buoy_*), and NWS warnings issued (nws_*, "
                "US zone, 2005 on). Thresholds: small craft 20 kn (Canada) / 23 kn (US), gale 34 kn, storm 48 kn.",
        "generated": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()), "areas": summary}, indent=1))
    log("done")


if __name__ == "__main__":
    main()
