#!/usr/bin/env python3
"""Candidate landing beaches where the route has no past campsite in reach.

Reads OpenStreetMap (Overpass) for beaches, shingle shores, campsites and slipways within NEAR_NM of the Stage 2 route, keeps
those more than FAR_NM from every past campsite (data/campsites.json, the stops past teams made), and writes
data/landings.json. For checking them by eye it also writes contact sheets of Esri satellite imagery around each candidate
(data/landings/sheets/*.jpg, 16 per sheet, ~600 m across each, numbered as in landings.json).
"""
from __future__ import annotations

import io
import json
import math
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent.parent
NEAR_NM, FAR_NM = 3.0, 3.0
BBOX = (48.3, -131.9, 55.6, -123.0)
ENDPOINTS = ["https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter",
             "https://maps.mail.ru/osm/tools/overpass/api/interpreter"]
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
B = ",".join(map(str, BBOX))
QUERY = f"""[out:json][timeout:300];
(nwr["natural"~"^(beach|shingle)$"]({B}); nwr["tourism"="camp_site"]({B}); nwr["leisure"="slipway"]({B});
 nwr["seamark:type"="small_craft_facility"]({B}););
out center tags;"""
IMG = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
Z = 17


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, data=None, timeout=60):
    with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=timeout) as r:
        return r.read()


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


def nm(a, b):
    return math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0])))


def tile_xy(lat, lon, z):
    n = 2 ** z
    x = (lon + 180) / 360 * n
    y = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
    return x, y


def chip(lat, lon, size=256):
    """~600 m of imagery centred on lat, lon at zoom Z, scaled to size px, with a crosshair."""
    fx, fy = tile_xy(lat, lon, Z)
    x0, y0 = int(fx) - 1, int(fy) - 1
    im = Image.new("RGB", (768, 768), (20, 30, 40))
    for dx in range(3):
        for dy in range(3):
            try:
                im.paste(Image.open(io.BytesIO(get(IMG.format(z=Z, x=x0 + dx, y=y0 + dy)))).convert("RGB"), (dx * 256, dy * 256))
            except Exception as e:
                log(f"  tile {e}")
    cx, cy = (fx - x0) * 256, (fy - y0) * 256
    half = 384
    im = im.crop((int(cx - half), int(cy - half), int(cx + half), int(cy + half))) if 0 <= cx - half and cx + half <= 768 else im
    im = im.resize((size, size))
    d = ImageDraw.Draw(im)
    c = size // 2
    d.line((c - 14, c, c - 5, c), fill=(255, 60, 60), width=2); d.line((c + 5, c, c + 14, c), fill=(255, 60, 60), width=2)
    d.line((c, c - 14, c, c - 5), fill=(255, 60, 60), width=2); d.line((c, c + 5, c, c + 14), fill=(255, 60, 60), width=2)
    return im


def main():
    route = json.loads((ROOT / "data" / "route.json").read_text())
    pts = [p for leg in route["legs"] if leg["key"] != "stage1" for p in leg["points"]]
    dense = []
    for a, b in zip(pts, pts[1:]):
        n = max(1, int(nm(a, b) / 0.5))
        dense += [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n)]
    camps = json.loads((ROOT / "data" / "campsites.json").read_text())
    els = overpass()
    log(f"{len(els)} candidate features in the box")
    out = []
    for e in els:
        la, lo = (e["lat"], e["lon"]) if "lat" in e else (e["center"]["lat"], e["center"]["lon"])
        dr = min(nm((la, lo), p) for p in dense)
        if dr > NEAR_NM or min(nm((la, lo), c) for c in camps) <= FAR_NM:
            continue
        t = e.get("tags", {})
        kind = t.get("natural") or ("campsite" if t.get("tourism") == "camp_site" else "slipway" if t.get("leisure") else "small craft")
        if any(nm((la, lo), (o["lat"], o["lon"])) < 0.3 for o in out):
            continue
        out.append({"lat": round(la, 5), "lon": round(lo, 5), "kind": kind, "name": t.get("name", ""), "osm": f"{e['type']}/{e['id']}",
                    "route_nm": round(dr, 2)})
    out.sort(key=lambda x: (x["lat"], x["lon"]))
    for i, o in enumerate(out):
        o["n"] = i + 1
    (ROOT / "data" / "landings.json").write_text(json.dumps({
        "source": "OpenStreetMap contributors (ODbL); imagery Esri World Imagery", "near_nm": NEAR_NM, "far_from_campsites_nm": FAR_NM,
        "generated": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()), "landings": out}, indent=0))
    log(f"{len(out)} candidates more than {FAR_NM} nm from a past campsite and within {NEAR_NM} nm of the route")
    sd = ROOT / "data" / "landings" / "sheets"
    sd.mkdir(parents=True, exist_ok=True)
    for f in sd.glob("*.jpg"):
        f.unlink()
    for s in range(0, len(out), 16):
        sheet = Image.new("RGB", (4 * 260, 4 * 278), (255, 255, 255))
        d = ImageDraw.Draw(sheet)
        for k, o in enumerate(out[s:s + 16]):
            x, y = (k % 4) * 260, (k // 4) * 278
            sheet.paste(chip(o["lat"], o["lon"]), (x + 2, y + 2))
            d.text((x + 4, y + 260), f"#{o['n']} {o['kind']} {o['lat']:.4f},{o['lon']:.4f}", fill=(0, 0, 0))
        sheet.save(sd / f"sheet{s // 16 + 1:02d}.jpg", quality=72)
        log(f"  sheet {s // 16 + 1}")


if __name__ == "__main__":
    main()
