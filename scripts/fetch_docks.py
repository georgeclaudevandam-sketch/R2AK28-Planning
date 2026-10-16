#!/usr/bin/env python3
"""Docks, floats and marinas near the Stage 2 route, where a boat can tie up for a night or a storm instead of landing on a beach.

Reads OpenStreetMap (Overpass) for marinas, piers and floating docks, harbours and small-craft facilities within NEAR_NM of
the Stage 2 route (south, inside and outside legs) and writes data/docks.json. Ferry terminals, industrial and fuel-only
piers are left out. Features within MERGE_NM of each other are merged into one site (a marina's piers are one place), keeping
the most useful kind (marina, then float, then harbour, then dock) and a name if any has one. For checking them by eye it also
writes contact sheets of Esri satellite imagery around each site (data/docks/sheets/*.jpg, 16 per sheet, ~400 m across).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_landings import chip, get, log, nm, ENDPOINTS  # noqa: E402

import urllib.parse  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
NEAR_NM, MERGE_NM = 2.0, 0.2
BBOX = (48.3, -131.9, 55.6, -123.0)
B = ",".join(map(str, BBOX))
QUERY = f"""[out:json][timeout:300];
(nwr["leisure"="marina"]({B}); nwr["man_made"="pier"]({B}); nwr["seamark:type"~"^(harbour|small_craft_facility)$"]({B});
 nwr["harbour"="yes"]({B}); nwr["mooring"]["mooring"!="no"]({B}););
out center tags;"""
RANK = {"marina": 0, "float": 1, "harbour": 2, "dock": 3}


def overpass():
    data = urllib.parse.urlencode({"data": QUERY}).encode()
    for i in range(6):
        url = ENDPOINTS[i % len(ENDPOINTS)]
        try:
            return json.loads(get(url, data, 400))["elements"]
        except Exception as e:
            log(f"  {url}: {type(e).__name__} {e}")
            time.sleep(15 * (i + 1))
    sys.exit("Overpass failed")


def kind_of(t):
    """None for places a small boat cannot use (ferries, industry, private, fuel only)."""
    if t.get("amenity") == "ferry_terminal" or "ferry" in t.get("name", "").lower() or t.get("public_transport"):
        return None
    if t.get("industrial") or t.get("landuse") == "industrial" or t.get("access") in ("private", "no"):
        return None
    if t.get("seamark:harbour:category") in ("ferry_terminal", "naval_base", "container_terminal", "bulk_terminal", "roro_terminal"):
        return None
    if t.get("leisure") == "marina":
        return "marina"
    if t.get("floating") == "yes" or "float" in t.get("name", "").lower() or t.get("seamark:type") == "small_craft_facility":
        return "float"
    if t.get("seamark:type") == "harbour" or t.get("harbour") == "yes":
        return "harbour"
    if t.get("man_made") == "pier" or t.get("mooring"):
        return "dock"
    return None


def main():
    route = json.loads((ROOT / "data" / "route.json").read_text())
    dense = []
    for leg in route["legs"]:
        if leg["key"] == "stage1":
            continue
        pts = leg["points"]
        for a, b in zip(pts, pts[1:]):
            n = max(1, int(nm(a, b) / 0.5))
            dense += [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n)]
    els = overpass()
    log(f"{len(els)} features in the box")
    out = []
    for e in els:
        la, lo = (e["lat"], e["lon"]) if "lat" in e else (e["center"]["lat"], e["center"]["lon"])
        t = e.get("tags", {})
        k = kind_of(t)
        if not k:
            continue
        dr = min(nm((la, lo), p) for p in dense)
        if dr > NEAR_NM:
            continue
        o = next((o for o in out if nm((la, lo), (o["lat"], o["lon"])) < MERGE_NM), None)
        if o:
            if RANK[k] < RANK[o["kind"]]:
                o["kind"] = k
            if not o["name"] and t.get("name"):
                o["name"] = t["name"]
            o["osm"].append(f"{e['type']}/{e['id']}")
            continue
        out.append({"lat": round(la, 5), "lon": round(lo, 5), "kind": k, "name": t.get("name", ""), "osm": [f"{e['type']}/{e['id']}"],
                    "route_nm": round(dr, 2)})
    out.sort(key=lambda x: (x["lat"], x["lon"]))
    for i, o in enumerate(out):
        o["n"] = i + 1
    (ROOT / "data" / "docks.json").write_text(json.dumps({
        "note": "Marinas, floats, harbours and docks within 2 nm of the Stage 2 route from OpenStreetMap (ODbL). Positions are the "
                "middle of the mapped feature. Whether a visiting rowboat may tie up is not known: check before relying on one. "
                "For planning, not navigation. Built by scripts/fetch_docks.py.",
        "near_nm": NEAR_NM, "generated": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()), "sites": out}, indent=0))
    log(f"{len(out)} sites within {NEAR_NM} nm of the route")
    sd = ROOT / "data" / "docks" / "sheets"
    sd.mkdir(parents=True, exist_ok=True)
    for f in sd.glob("*.jpg"):
        f.unlink()
    for s in range(0, len(out), 16):
        sheet = Image.new("RGB", (4 * 260, 4 * 278), (255, 255, 255))
        d = ImageDraw.Draw(sheet)
        for k, o in enumerate(out[s:s + 16]):
            x, y = (k % 4) * 260, (k // 4) * 278
            sheet.paste(chip(o["lat"], o["lon"]), (x + 2, y + 2))
            d.text((x + 4, y + 260), f"#{o['n']} {o['kind']} {o['name'][:22]}", fill=(0, 0, 0))
        sheet.save(sd / f"sheet{s // 16 + 1:02d}.jpg", quality=72)
        log(f"  sheet {s // 16 + 1}")


if __name__ == "__main__":
    main()
