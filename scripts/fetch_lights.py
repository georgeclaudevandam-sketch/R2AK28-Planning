#!/usr/bin/env python3
"""Fetch navigation lights near the route from OpenStreetMap (OpenSeaMap seamark tags) via the Overpass API.

Writes data/lights.json: every light within NEAR_NM of a route leg, with its characteristic per sector (character, colours,
period, group and, where mapped, the exact on/off sequence) so the map can flash it on its real rhythm.
"""
from __future__ import annotations

import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NEAR_NM = 3.0
BBOX = (47.9, -131.9, 55.6, -122.5)   # south, west, north, east
ENDPOINTS = ["https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter",
             "https://maps.mail.ru/osm/tools/overpass/api/interpreter"]
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
QUERY = f"""[out:json][timeout:300];
(node[~"^seamark:light(:[0-9]+)?:(character|colour)$"~"."]({BBOX[0]},{BBOX[1]},{BBOX[2]},{BBOX[3]}););
out body;"""


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def overpass():
    data = urllib.parse.urlencode({"data": QUERY}).encode()
    for i in range(6):
        url = ENDPOINTS[i % len(ENDPOINTS)]
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=400) as r:
                return json.load(r)["elements"]
        except Exception as e:
            log(f"  {url}: {type(e).__name__} {e}")
            time.sleep(15 * (i + 1))
    sys.exit("Overpass failed")


def num(v):
    try:
        return float(str(v).strip().rstrip("s"))
    except (TypeError, ValueError):
        return None


def sectors(tags):
    """One entry per light sector: seamark:light:N:* (or the unnumbered seamark:light:* for a single light)."""
    idx = sorted({int(m.group(1)) for k in tags if (m := re.match(r"seamark:light:(\d+):", k))})
    pre = [f"seamark:light:{i}:" for i in idx] or ["seamark:light:"]
    out = []
    for p in pre:
        g = lambda k: tags.get(p + k) or (tags.get("seamark:light:" + k) if p != "seamark:light:" else None)
        ch, col = g("character"), g("colour")
        if not ch and not col:
            continue
        s = {"ch": ch or "F", "col": (col or "white").split(";")[0]}
        for k, key in (("period", "per"), ("range", "rng"), ("sector_start", "s0"), ("sector_end", "s1")):
            v = num(g(k))
            if v is not None:
                s[key] = v
        for k, key in (("group", "grp"), ("sequence", "seq")):
            if g(k):
                s[key] = g(k)
        out.append(s)
    return out


def main():
    route = json.loads((ROOT / "data" / "route.json").read_text())
    pts = [p for leg in route["legs"] for p in leg["points"]]
    # dense route points for the distance test
    dense = []
    for a, b in zip(pts, pts[1:]):
        n = max(1, int(math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0]))) / 0.5))
        dense += [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n)]
    cells = {}
    for la, lo in dense:
        cells.setdefault((int(la * 10), int(lo * 10)), []).append((la, lo))

    def near(la, lo):
        k = math.cos(math.radians(la))
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                for a, b in cells.get((int(la * 10) + dy, int(lo * 10) + dx), []):
                    if math.hypot((a - la) * 60, (b - lo) * 60 * k) <= NEAR_NM:
                        return True
        return False

    els = overpass()
    log(f"{len(els)} lit seamarks in the box")
    lights = []
    for e in els:
        t = e.get("tags", {})
        if not near(e["lat"], e["lon"]):
            continue
        sec = sectors(t)
        if not sec:
            continue
        lights.append({"id": e["id"], "lat": round(e["lat"], 5), "lon": round(e["lon"], 5),
                       "name": t.get("seamark:name") or t.get("name") or "", "type": t.get("seamark:type", ""), "s": sec})
    lights.sort(key=lambda x: (-x["lat"], x["lon"]))
    (ROOT / "data" / "lights.json").write_text(json.dumps({
        "source": "OpenStreetMap contributors (OpenSeaMap seamark tags), ODbL", "near_nm": NEAR_NM,
        "generated": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()), "lights": lights}, separators=(",", ":")))
    seq = sum(1 for l in lights if any("seq" in s for s in l["s"]))
    log(f"wrote {len(lights)} lights within {NEAR_NM} nm of the route ({seq} with an exact sequence)")


if __name__ == "__main__":
    main()
