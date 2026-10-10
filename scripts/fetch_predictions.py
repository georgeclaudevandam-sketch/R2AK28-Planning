#!/usr/bin/env python3
"""Fetch tide and current predictions for the R2AK 2028 planning map.

Sources
  NOAA CO-OPS  https://api.tidesandcurrents.noaa.gov  (US stations, and NOAA's
               subordinate current stations in northern BC)
  CHS IWLS     https://api-iwls.dfo-mpo.gc.ca          (Canadian stations)

CHS publishes current predictions about a year ahead. When a CHS current
station has no official events for the 2028 window yet, this script derives
them: it fits harmonic constituents (UTide) to CHS's own published 2027
prediction series for that station, checks the fit against CHS's Oct-Dec 2027
events, and then predicts the 2028 window. Those stations are flagged
"derived" in the output, with the validation numbers, and should be replaced by
official CHS tables once they are published (re-run this script).

Output: data/predictions.json
"""
from __future__ import annotations

import json
import math
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STATIONS = ROOT / "data" / "stations.json"
OUT = ROOT / "data" / "predictions.json"

NOAA = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
NOAA_MD = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations"
IWLS = "https://api-iwls.dfo-mpo.gc.ca/api/v1"
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}

# Year of published CHS predictions used to derive constituents when the
# target year is not published yet.
FIT_YEAR = 2027


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get_json(url: str, params: dict | None = None, tries: int = 4):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            body = e.read().decode(errors="replace")[:300]
            if e.code in (400, 404):
                raise RuntimeError(f"HTTP {e.code} for {url}: {body}") from None
            log(f"  HTTP {e.code}, retry {i + 1}: {url}")
        except Exception as e:  # network hiccup
            log(f"  {type(e).__name__} {e}, retry {i + 1}: {url}")
        time.sleep(2 ** (i + 1))
    raise RuntimeError(f"giving up on {url}")


def iso(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")


def parse_utc(s: str) -> datetime:
    s = s.replace("Z", "").replace("T", " ")
    return datetime.strptime(s[:16], "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)


def chunks(start: datetime, end: datetime, days: int):
    t = start
    while t < end:
        u = min(end, t + timedelta(days=days))
        yield t, u
        t = u


def label_hilo(events: list[dict]) -> list[dict]:
    """CHS hilo values carry no H/L flag: label by comparing neighbours."""
    for i, e in enumerate(events):
        prev = events[i - 1]["h"] if i > 0 else None
        nxt = events[i + 1]["h"] if i + 1 < len(events) else None
        ref = [x for x in (prev, nxt) if x is not None]
        e["type"] = "H" if all(e["h"] > x for x in ref) else "L"
    return events


# ---------------------------------------------------------------- NOAA

def noaa_tide(st, start, end):
    out = []
    for a, b in chunks(start, end, 31):
        d = get_json(NOAA, {
            "product": "predictions", "application": "r2ak28-planning",
            "begin_date": a.strftime("%Y%m%d %H:%M"), "end_date": b.strftime("%Y%m%d %H:%M"),
            "datum": "MLLW", "station": st["id"], "time_zone": "gmt", "units": "metric",
            "interval": "hilo", "format": "json"})
        if "error" in d:
            raise RuntimeError(d["error"])
        for p in d["predictions"]:
            out.append({"t": iso(parse_utc(p["t"])), "type": p["type"], "h": round(float(p["v"]), 2)})
    return {"datum": "MLLW", "units": "m", "events": dedupe(out)}


def noaa_current(st, start, end):
    out, meta = [], {}
    params = {
        "product": "currents_predictions", "application": "r2ak28-planning",
        "station": st["id"], "time_zone": "gmt", "units": "english",
        "interval": "MAX_SLACK", "format": "json"}
    if st.get("bin"):
        params["bin"] = st["bin"]
    for a, b in chunks(start, end, 31):
        p = dict(params, begin_date=a.strftime("%Y%m%d %H:%M"), end_date=b.strftime("%Y%m%d %H:%M"))
        d = get_json(NOAA, p)
        if "error" in d and "bin" not in params:
            # multi-bin station: ask the metadata API for its bins and use the first
            bins = get_json(f"{NOAA_MD}/{st['id']}/bins.json", {"units": "english"})
            params["bin"] = str(bins["bins"][0]["num"])
            log(f"  {st['id']}: using bin {params['bin']}")
            d = get_json(NOAA, dict(params, begin_date=p["begin_date"], end_date=p["end_date"]))
        if "error" in d:
            raise RuntimeError(d["error"])
        for c in d["current_predictions"]["cp"]:
            typ = c.get("Type", "").lower()
            v = float(c["Velocity_Major"])
            meta = {"flood_dir": c.get("meanFloodDir"), "ebb_dir": c.get("meanEbbDir")}
            out.append({"t": iso(parse_utc(c["Time"])), "type": typ, "kn": round(abs(v), 2)})
    return {"units": "kn", **meta, "bin": params.get("bin"), "events": dedupe(out)}


# ---------------------------------------------------------------- CHS

_chs_ids: dict[str, dict] = {}


def chs_station(code: str) -> dict:
    if code not in _chs_ids:
        d = get_json(f"{IWLS}/stations", {"code": code})
        if not d:
            raise RuntimeError(f"CHS station code {code} not found")
        _chs_ids[code] = d[0]
    return _chs_ids[code]


def chs_series(code, series, start, end, days=7, resolution=None):
    sid = chs_station(code)["id"]
    out = []
    for a, b in chunks(start, end, days):
        p = {"time-series-code": series, "from": a.strftime("%Y-%m-%dT%H:%M:%SZ"),
             "to": b.strftime("%Y-%m-%dT%H:%M:%SZ")}
        if resolution:
            p["resolution"] = resolution
        out += get_json(f"{IWLS}/stations/{sid}/data", p)
        time.sleep(0.3)
    return out


def chs_tide(st, start, end):
    raw = chs_series(st["id"], "wlp-hilo", start, end, days=31)
    ev = [{"t": iso(parse_utc(r["eventDate"])), "h": round(float(r["value"]), 2)} for r in raw]
    ev = label_hilo(dedupe(ev))
    if not ev:
        raise RuntimeError("no CHS hilo predictions returned")
    return {"datum": "chart datum (LLWLT)", "units": "m", "events": ev}


def chs_events_to_list(raw):
    q = {"SLACK": "slack", "EXTREMA_FLOOD": "flood", "EXTREMA_EBB": "ebb"}
    ev = []
    for r in raw:
        typ = q.get(r.get("qualifier", ""), None)
        if typ is None:
            continue
        ev.append({"t": iso(parse_utc(r["eventDate"])), "type": typ, "kn": round(abs(float(r["value"])), 2)})
    return dedupe(ev)


def chs_current(st, start, end):
    raw = chs_series(st["id"], "wcp1-events", start, end, days=31)
    ev = chs_events_to_list(raw)
    span = (end - start).total_seconds() / 86400
    if len(ev) >= 2 * span:  # roughly complete: official tables exist
        return {"units": "kn", "events": ev, "method": "official"}
    log(f"  {st['name']}: {len(ev)} official events in window, deriving from {FIT_YEAR} harmonics")
    return chs_current_derived(st, start, end)


def chs_current_derived(st, start, end):
    import numpy as np
    import utide

    y0 = datetime(FIT_YEAR, 1, 1, tzinfo=timezone.utc)
    y1 = datetime(FIT_YEAR + 1, 1, 1, tzinfo=timezone.utc)
    # 15-minute speeds; try a monthly request first, fall back to weekly
    try:
        speed = chs_series(st["id"], "wcsp1", y0, y1, days=31, resolution="FIFTEEN_MINUTES")
    except RuntimeError as e:
        log(f"  monthly 15-min request refused ({e}); using weekly requests")
        speed = chs_series(st["id"], "wcsp1", y0, y1, days=7, resolution="FIFTEEN_MINUTES")
    events = chs_events_to_list(chs_series(st["id"], "wcp1-events", y0, y1, days=31))
    if len(speed) < 1000 or len(events) < 100:
        raise RuntimeError(f"not enough {FIT_YEAR} CHS data to derive constituents")

    t = np.array([parse_utc(r["eventDate"]).timestamp() for r in speed])
    s = np.array([abs(float(r["value"])) for r in speed])
    order = np.argsort(t)
    t, s = t[order], s[order]
    # sign the unsigned speed: flood +, ebb -, using the extreme between each pair of slacks
    et = np.array([parse_utc(e["t"]).timestamp() for e in events])
    sign = np.zeros_like(s)
    slack_t = [et[i] for i, e in enumerate(events) if e["type"] == "slack"]
    ext = [(et[i], 1.0 if e["type"] == "flood" else -1.0) for i, e in enumerate(events) if e["type"] != "slack"]
    for a, b in zip(slack_t, slack_t[1:]):
        inside = [sg for tt, sg in ext if a < tt < b]
        if inside:
            sign[(t >= a) & (t < b)] = inside[0]
    keep = sign != 0
    t, u = t[keep], (s * sign)[keep]
    lat = st["lat"]
    as_dt64 = lambda secs: (np.asarray(secs) * 1e9).astype("int64").astype("datetime64[ns]")

    def fit(mask):
        return utide.solve(as_dt64(t[mask]), u[mask], lat=lat, method="ols", conf_int="none",
                           trend=False, nodal=True, verbose=False)

    def predict(coef, a: datetime, b: datetime):
        tt = np.arange(a.timestamp(), b.timestamp(), 60.0)
        return tt, utide.reconstruct(as_dt64(tt), coef, verbose=False).h

    # validation: fit Jan-Sep, compare Oct-Dec with CHS's own events
    split = datetime(FIT_YEAR, 10, 1, tzinfo=timezone.utc)
    coef_v = fit(t < split.timestamp())
    vt, vu = predict(coef_v, split, y1)
    derived_v = extract_events(vt, vu)
    official_v = [e for e in events if parse_utc(e["t"]) >= split]
    validation = compare_events(derived_v, official_v)
    log(f"  validation {st['name']}: {validation}")

    coef = fit(np.ones_like(t, dtype=bool))
    pt, pu = predict(coef, start, end)
    ev = extract_events(pt, pu)
    return {"units": "kn", "events": ev, "method": "derived",
            "derived_from": f"UTide fit to CHS {FIT_YEAR} published predictions (15-min speeds, signed by CHS events)",
            "constituents": int(len(coef["name"])), "validation": validation}


def extract_events(tt, uu):
    """Slack = zero crossing; flood/ebb maximum = largest |u| between crossings."""
    import numpy as np

    ev = []
    z = np.where(np.sign(uu[:-1]) != np.sign(uu[1:]))[0]
    for a, b in zip(z, z[1:]):
        seg = uu[a + 1:b + 1]
        if len(seg) < 30:  # ignore crossings less than 30 min apart (noise at weak slacks)
            continue
        k = int(np.argmax(np.abs(seg))) + a + 1
        ev.append({"t": iso(datetime.fromtimestamp(tt[k], timezone.utc)),
                   "type": "flood" if uu[k] > 0 else "ebb", "kn": round(float(abs(uu[k])), 2)})
    for i in z:
        f = uu[i] / (uu[i] - uu[i + 1])
        ts = tt[i] + f * 60.0
        ev.append({"t": iso(datetime.fromtimestamp(ts, timezone.utc)), "type": "slack", "kn": 0.0})
    return dedupe(ev)


def compare_events(derived, official):
    import numpy as np

    def match(kind):
        d = [(parse_utc(e["t"]).timestamp(), e["kn"]) for e in derived if (e["type"] == "slack") == (kind == "slack")]
        if not d:
            return [], []
        dt = np.array([x[0] for x in d])
        terr, verr = [], []
        for e in official:
            if (e["type"] == "slack") != (kind == "slack"):
                continue
            ot = parse_utc(e["t"]).timestamp()
            i = int(np.argmin(np.abs(dt - ot)))
            if abs(dt[i] - ot) < 3 * 3600:
                terr.append(abs(dt[i] - ot) / 60.0)
                if kind != "slack":
                    verr.append(abs(d[i][1] - e["kn"]))
        return terr, verr

    st, _ = match("slack")
    mt, mv = match("max")
    q = lambda a, p: round(float(np.percentile(a, p)), 1) if len(a) else None
    return {"period": "Oct-Dec %d (held out)" % FIT_YEAR,
            "slack_minutes_median": q(st, 50), "slack_minutes_p90": q(st, 90),
            "max_minutes_median": q(mt, 50), "max_knots_median": q(mv, 50), "max_knots_p90": q(mv, 90),
            "matched_slacks": len(st), "official_events": len(official)}


def dedupe(ev):
    seen, out = set(), []
    for e in sorted(ev, key=lambda e: e["t"]):
        k = (e["t"], e.get("type"))
        if k not in seen:
            seen.add(k)
            out.append(e)
    return out


# ---------------------------------------------------------------- main

def main():
    cfg = json.loads(STATIONS.read_text())
    start = datetime.fromisoformat(cfg["window"]["start"].replace("Z", "+00:00"))
    end = datetime.fromisoformat(cfg["window"]["end"].replace("Z", "+00:00"))
    only = set(sys.argv[1:])
    prev = json.loads(OUT.read_text())["stations"] if OUT.exists() else {}
    result, failures = {}, []
    for st in cfg["stations"]:
        if only and st["key"] not in only:
            if st["key"] in prev:
                result[st["key"]] = prev[st["key"]]
            continue
        log(f"{st['source'].upper()} {st['kind']} {st['id']} {st['name']}")
        fn = {("noaa", "tide"): noaa_tide, ("noaa", "current"): noaa_current,
              ("chs", "tide"): chs_tide, ("chs", "current"): chs_current}[(st["source"], st["kind"])]
        try:
            data = fn(st, start, end)
            data.setdefault("method", "official")
            result[st["key"]] = data
            log(f"  {len(data['events'])} events ({data['method']})")
        except Exception as e:  # keep going; report at the end
            log(f"  FAILED: {e}")
            failures.append(f"{st['key']}: {e}")
            if st["key"] in prev:
                result[st["key"]] = prev[st["key"]]
    OUT.write_text(json.dumps({
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "window": cfg["window"], "failures": failures, "stations": result,
    }, separators=(",", ":")))
    log(f"wrote {OUT} ({len(result)} stations, {len(failures)} failures)")
    for f in failures:
        log("  " + f)


if __name__ == "__main__":
    main()
