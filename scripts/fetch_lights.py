#!/usr/bin/env python3
"""Fetch navigation lights near the route from the official electronic charts (CHS and NOAA ENC chart services).

Asks each Maritime Chart Service to identify the LIGHTS features (S-57) in boxes covering the route corridor, and writes
data/lights.json: every light within NEAR_NM of a route leg, with its characteristic per sector (character, colour, period,
group and sector bearings) so the map can flash it on its real rhythm.
"""
from __future__ import annotations

import json
import math
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NEAR_NM = 5.0
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
SVC = {
    "chs": ("https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/chs/ENC_MaritimeChartService/MapServer/exts/MaritimeChartService/MapServer",
            lambda la, lo: 48.15 < la < 54.8),
    "noaa": ("https://gis.charttools.noaa.gov/arcgis/rest/services/MCS/NOAAChartDisplay/MapServer/exts/MaritimeChartService/MapServer",
             lambda la, lo: la < 49.05 or la > 54.55),
}
DLA, DLO = 0.06, 0.1   # box size, degrees
# S-57 light characteristic (LITCHR) as words -> IALA abbreviation; longest first so "very quick" wins over "quick"
CHR = [("ultra quick", "UQ"), ("interrupted very quick", "IVQ"), ("interrupted quick", "IQ"), ("very quick", "VQ"), ("quick", "Q"),
       ("long-flash", "LFl"), ("long flash", "LFl"), ("fixed and flash", "FFl"), ("flash", "Fl"), ("isophase", "Iso"), ("occult", "Oc"),
       ("morse", "Mo"), ("fixed", "F")]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return json.load(r)
        except Exception as e:
            log(f"  {type(e).__name__} {e}")
            time.sleep(5 * (i + 1))
    return None


def nm(a, b):
    return math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0])))


def attr(a, *names):
    low = {k.strip().lower(): v for k, v in a.items()}
    for n in names:
        v = low.get(n.lower())
        if v not in (None, "", "Unknown"):
            return str(v)
    return None


def num(v):
    try:
        return float(str(v).split()[0])
    except (TypeError, ValueError, IndexError):
        return None


def light(a):
    raw = (attr(a, "LITCHR", "Light characteristic") or "fixed").lower()
    ch = next((ab for w, ab in CHR if w in raw), "F")
    if "alternat" in raw:
        ch = "Al." + ch
    if "long-flash" in raw and ch in ("Q", "VQ"):
        ch += "+LFl"
    s = {"ch": ch, "col": (attr(a, "COLOUR", "Colour") or "white").replace(",", ";")}
    grp = (attr(a, "SIGGRP", "Signal group") or "").strip("()")
    if grp and grp not in ("1",) and any(c.isdigit() for c in grp):
        s["grp"] = grp
    for key, names in (("per", ("SIGPER", "Signal period")), ("s0", ("SECTR1", "Sector limit one")), ("s1", ("SECTR2", "Sector limit two")),
                       ("rng", ("VALNMR", "Value of nominal range"))):
        v = num(attr(a, *names))
        if v is not None:
            s[key] = v
    sq = attr(a, "SIGSEQ", "Signal sequence")
    if sq:
        s["seq"] = sq
    return s


def main():
    route = json.loads((ROOT / "data" / "route.json").read_text())
    dense = []
    for leg in route["legs"]:
        p = leg["points"]
        for a, b in zip(p, p[1:]):
            n = max(1, int(nm(a, b) / 0.5))
            dense += [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n + 1)]
    pad = NEAR_NM / 60
    boxes = sorted({(math.floor(la / DLA), math.floor(lo / DLO)) for la, lo in dense
                    for la2 in (la - pad, la, la + pad) for lo2 in (lo - pad * 1.6, lo, lo + pad * 1.6)
                    for la, lo in [(la2, lo2)]})
    log(f"{len(boxes)} boxes")
    jobs = [(k, i, j) for i, j in boxes for k, (_, inside) in SVC.items() if inside((i + 0.5) * DLA, (j + 0.5) * DLO)]

    def ident(job):
        k, i, j = job
        x0, y0, x1, y1 = j * DLO, i * DLA, (j + 1) * DLO, (i + 1) * DLA
        q = {"geometry": json.dumps({"xmin": x0, "ymin": y0, "xmax": x1, "ymax": y1, "spatialReference": {"wkid": 4326}}),
             "geometryType": "esriGeometryEnvelope", "sr": "4326", "layers": "all:6", "tolerance": "1", "mapExtent": f"{x0},{y0},{x1},{y1}",
             "imageDisplay": "1200,900,96", "returnGeometry": "true", "f": "json"}
        r = get(SVC[k][0] + "/identify?" + urllib.parse.urlencode(q)) or {}
        return k, [x for x in r.get("results", []) if (x.get("attributes") or {}).get("objectType") == "LIGHTS"]

    feats, seen = {}, set()
    with ThreadPoolExecutor(6) as ex:
        for n, (k, res) in enumerate(ex.map(ident, jobs)):
            for x in res:
                a, g = x["attributes"], x.get("geometry") or {}
                if "x" not in g:
                    continue
                key = attr(a, "LNAM") or f"{k}{a.get('rcid')}{g['x']:.5f}{g['y']:.5f}"
                if key in seen:
                    continue
                seen.add(key)
                if "flood" in (attr(a, "CATLIT", "Category of light") or "").lower():
                    continue
                ll = (round(g["y"], 5), round(g["x"], 5))
                f = feats.setdefault(ll, {"lat": ll[0], "lon": ll[1], "name": attr(a, "OBJNAM", "Object name") or "", "src": k, "s": []})
                f["s"].append(light(a))
            if n % 50 == 0:
                log(f"  {n}/{len(jobs)} requests, {len(feats)} lights")
    lights = [f for f in feats.values() if min(nm((f["lat"], f["lon"]), p) for p in dense[::2]) <= NEAR_NM]
    for f in lights:   # sectors in bearing order, one entry per distinct sector
        f["s"] = sorted({json.dumps(s, sort_keys=True) for s in f["s"]}, key=lambda t: json.loads(t).get("s0", -1))
        f["s"] = [json.loads(t) for t in f["s"]]
    lights.sort(key=lambda x: (-x["lat"], x["lon"]))
    (ROOT / "data" / "lights.json").write_text(json.dumps({
        "source": "CHS and NOAA electronic navigational charts (S-57 LIGHTS) via their Maritime Chart Services. Not for navigation.",
        "near_nm": NEAR_NM, "generated": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()), "lights": lights}, separators=(",", ":")))
    log(f"wrote {len(lights)} lights within {NEAR_NM} nm of the route ({sum(f['src'] == 'chs' for f in lights)} CHS, "
        f"{sum(f['src'] == 'noaa' for f in lights)} NOAA)")


if __name__ == "__main__":
    main()
