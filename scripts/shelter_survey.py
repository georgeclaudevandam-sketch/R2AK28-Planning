#!/usr/bin/env python3
"""Storm-shelter survey of one stretch of the route (default: Tolmie Channel and Finlayson Channel north of Klemtu, where the
storm page finds no shelter in reach).

Gathers, within NEAR_NM of the route between latitudes LAT0 and LAT1:
  - chart features from the CHS ENC (Maritime Chart Service identify): anchorages, moorings, small-craft and harbour
    facilities, piers and other shoreline constructions, buildings, landmarks, seabed nature, drying areas;
  - OpenStreetMap features: creek and river mouths (deltas often give a gravel landing), beaches, bays, piers, buildings,
    places, seamarks;
  - Esri satellite chips (~400 m) of the shoreline every ~0.2 nm and of every feature above, on numbered contact sheets.

Writes data/shelter_survey/features.json, data/shelter_survey/index.json and data/shelter_survey/sheets/*.jpg.
Usage: python scripts/shelter_survey.py [LAT0 LAT1]
"""
from __future__ import annotations

import io
import json
import math
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "shelter_survey"
LAT0, LAT1 = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (52.54, 52.93)
NEAR_NM = 3.0
COAST_STEP_NM = 0.2
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
CHS = "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/chs/ENC_MaritimeChartService/MapServer/exts/MaritimeChartService/MapServer"
OVERPASS = ["https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter",
            "https://maps.mail.ru/osm/tools/overpass/api/interpreter"]
IMG = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
# S-57 object classes worth a look for shelter (others are counted only)
KEEP = {"ACHARE", "ACHBRT", "MORFAC", "SMCFAC", "HRBFAC", "HRBARE", "SLCONS", "BUISGL", "BUAARE", "LNDMRK", "SBDARE", "PILPNT",
        "CRANES", "OFSPLF", "PONTON", "HULKES", "FSHFAC", "CGUSTA", "RSCSTA", "LNDRGN", "SEAARE", "RIVERS", "CANALS", "LAKARE",
        "DEPARE", "CTNARE", "MARCUL", "RUNWAY", "AIRARE", "ROADWY", "CHKPNT", "SILTNK", "LIGHTS", "BCNSPP", "DAYMAR"}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def fetch(url, data=None, timeout=120, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=timeout) as r:
                return r.read()
        except Exception as e:
            log(f"  {type(e).__name__} {e}")
            time.sleep(5 * (i + 1))
    return None


def nm(a, b):
    return math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0])))


def route_dense():
    route = json.loads((ROOT / "data" / "route.json").read_text())
    dense = []
    for leg in route["legs"]:
        if leg["key"] in ("stage1", "outside"):
            continue
        p = leg["points"]
        for a, b in zip(p, p[1:]):
            n = max(1, int(nm(a, b) / 0.25))
            dense += [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n + 1)]
    return [q for q in dense if LAT0 - 0.06 <= q[0] <= LAT1 + 0.06]


def dist_route(ll, dense):
    return min(nm(ll, q) for q in dense)


# ---------- charts
def chart_features(bbox, dense):
    s, w, n, e = bbox
    dla, dlo = 0.03, 0.05
    boxes = [(s + i * dla, w + j * dlo) for i in range(math.ceil((n - s) / dla)) for j in range(math.ceil((e - w) / dlo))]
    boxes = [b for b in boxes if dist_route((b[0] + dla / 2, b[1] + dlo / 2), dense) <= NEAR_NM + 2]
    log(f"chart: {len(boxes)} boxes")

    def ident(b):
        y0, x0 = b
        y1, x1 = y0 + dla, x0 + dlo
        q = {"geometry": json.dumps({"xmin": x0, "ymin": y0, "xmax": x1, "ymax": y1, "spatialReference": {"wkid": 4326}}),
             "geometryType": "esriGeometryEnvelope", "sr": "4326", "layers": "all", "tolerance": "1", "mapExtent": f"{x0},{y0},{x1},{y1}",
             "imageDisplay": "1200,900,96", "returnGeometry": "true", "f": "json"}
        r = fetch(CHS + "/identify?" + urllib.parse.urlencode(q))
        try:
            return json.loads(r).get("results", []) if r else []
        except ValueError:
            return []

    counts, out, seen = {}, [], set()
    with ThreadPoolExecutor(6) as ex:
        for res in ex.map(ident, boxes):
            for x in res:
                a = x.get("attributes") or {}
                t = a.get("objectType") or x.get("layerName") or "?"
                counts[t] = counts.get(t, 0) + 1
                if t not in KEEP:
                    continue
                g = x.get("geometry") or {}
                if "x" in g:
                    ll = (g["y"], g["x"])
                elif g.get("rings") or g.get("paths"):
                    pts = [p for r in (g.get("rings") or g.get("paths")) for p in r]
                    ll = (sum(p[1] for p in pts) / len(pts), sum(p[0] for p in pts) / len(pts))
                else:
                    continue
                key = a.get("LNAM") or f"{t}{ll[0]:.5f}{ll[1]:.5f}"
                if key in seen:
                    continue
                seen.add(key)
                if t == "DEPARE":   # only drying areas (intertidal flats): DRVAL1 below chart datum
                    try:
                        if float(str(a.get("DRVAL1", a.get("Depth range value 1", "0"))).split()[0]) >= 0:
                            continue
                    except ValueError:
                        continue
                d = dist_route(ll, dense)
                if d > NEAR_NM or not (LAT0 <= ll[0] <= LAT1):
                    continue
                attrs = {k: v for k, v in a.items() if v not in (None, "", "Unknown", "Null") and k not in ("objectType",) and len(str(v)) < 200}
                out.append({"src": "chart", "type": t, "lat": round(ll[0], 5), "lon": round(ll[1], 5), "route_nm": round(d, 2), "attrs": attrs})
    log("chart object types: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])))
    return out, counts


# ---------- OpenStreetMap
def osm(bbox):
    s, w, n, e = bbox
    B = f"{s},{w},{n},{e}"
    q = f"""[out:json][timeout:300];
(way["natural"="coastline"]({B});
 way["waterway"~"^(stream|river|brook)$"]({B});
 nwr["natural"~"^(beach|shingle|bay|sand|wetland|spring)$"]({B});
 nwr["man_made"~"^(pier|jetty|breakwater|lighthouse|tower|mast)$"]({B});
 nwr["building"]({B}); nwr["place"]({B}); nwr["amenity"]({B}); nwr["leisure"]({B}); nwr["tourism"]({B});
 nwr["harbour"]({B}); nwr[~"^seamark:type$"~"."]({B}); nwr["landuse"]({B}););
out geom tags;"""
    data = urllib.parse.urlencode({"data": q}).encode()
    for i in range(6):
        r = fetch(OVERPASS[i % len(OVERPASS)], data, 400, tries=1)
        if r:
            try:
                return json.loads(r)["elements"]
            except ValueError:
                pass
        time.sleep(15 * (i + 1))
    log("Overpass failed")
    return []


def centre(e):
    if "lat" in e:
        return e["lat"], e["lon"]
    g = e.get("geometry") or []
    if g:
        return sum(p["lat"] for p in g) / len(g), sum(p["lon"] for p in g) / len(g)
    b = e.get("bounds")
    return ((b["minlat"] + b["maxlat"]) / 2, (b["minlon"] + b["maxlon"]) / 2) if b else None


# ---------- imagery
def tile_xy(lat, lon, z):
    n = 2 ** z
    return (lon + 180) / 360 * n, (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n


TILES = {}


def tile(z, x, y):
    k = (z, x, y)
    if k not in TILES:
        b = fetch(IMG.format(z=z, x=x, y=y), timeout=60, tries=3)
        TILES[k] = Image.open(io.BytesIO(b)).convert("RGB") if b else Image.new("RGB", (256, 256), (20, 30, 40))
    return TILES[k]


def chip(lat, lon, size=240, z=17):
    fx, fy = tile_xy(lat, lon, z)
    x0, y0 = int(fx) - 1, int(fy) - 1
    im = Image.new("RGB", (768, 768))
    for dx in range(3):
        for dy in range(3):
            im.paste(tile(z, x0 + dx, y0 + dy), (dx * 256, dy * 256))
    cx, cy = (fx - x0) * 256, (fy - y0) * 256
    im = im.crop((int(cx - 256), int(cy - 256), int(cx + 256), int(cy + 256))).resize((size, size))
    d = ImageDraw.Draw(im)
    c = size // 2
    for a, b in (((c - 14, c), (c - 5, c)), ((c + 5, c), (c + 14, c)), ((c, c - 14), (c, c - 5)), ((c, c + 5), (c, c + 14))):
        d.line((*a, *b), fill=(255, 60, 60), width=2)
    return im


def main():
    dense = route_dense()
    pad = (NEAR_NM + 0.5) / 60
    bbox = (LAT0 - pad, min(q[1] for q in dense) - pad * 1.7, LAT1 + pad, max(q[1] for q in dense) + pad * 1.7)
    log(f"route {len(dense)} points, box {bbox}")
    OUT.mkdir(parents=True, exist_ok=True)

    chart, counts = chart_features(bbox, dense)
    els = osm(bbox)
    log(f"OSM: {len(els)} elements")
    feats = list(chart)
    coast = []
    for e in els:
        t = e.get("tags", {})
        if t.get("natural") == "coastline":
            g = e.get("geometry") or []
            for a, b in zip(g, g[1:]):
                pa, pb = (a["lat"], a["lon"]), (b["lat"], b["lon"])
                n = max(1, int(nm(pa, pb) / 0.05))
                coast += [(pa[0] + (pb[0] - pa[0]) * k / n, pa[1] + (pb[1] - pa[1]) * k / n) for k in range(n)]
    log(f"coastline: {len(coast)} points")
    for e in els:
        t = e.get("tags", {})
        if t.get("natural") == "coastline":
            continue
        if t.get("waterway"):   # the mouth: the way's end nearest the shore
            g = e.get("geometry") or []
            if not g or not coast:
                continue
            ends = [(g[0]["lat"], g[0]["lon"]), (g[-1]["lat"], g[-1]["lon"])]
            m = min(ends, key=lambda p: min(nm(p, c) for c in coast[::3]))
            if min(nm(m, c) for c in coast) > 0.15:
                continue
            ll, kind = m, f"{t['waterway']} mouth"
        else:
            ll = centre(e)
            if not ll:
                continue
            kind = next((f"{k}={t[k]}" for k in ("natural", "man_made", "harbour", "seamark:type", "amenity", "leisure", "tourism", "place", "landuse") if k in t), "building" if "building" in t else "?")
        if not (LAT0 <= ll[0] <= LAT1):
            continue
        d = dist_route(ll, dense)
        if d > NEAR_NM:
            continue
        feats.append({"src": "osm", "type": kind, "name": t.get("name", ""), "lat": round(ll[0], 5), "lon": round(ll[1], 5), "route_nm": round(d, 2),
                      "osm": f"{e['type']}/{e['id']}", "tags": {k: v for k, v in t.items() if len(v) < 120}})
    feats.sort(key=lambda f: (-f["lat"], f["lon"]))
    for i, f in enumerate(feats):
        f["n"] = i + 1
    log(f"{len(feats)} features ({len(chart)} chart, {len(feats) - len(chart)} OSM)")

    # shoreline points every COAST_STEP_NM within NEAR_NM of the route, then a chip of each, then a chip of each feature
    shore = []
    for c in coast:
        if LAT0 <= c[0] <= LAT1 and all(nm(c, s) >= COAST_STEP_NM for s in shore[-40:]) and dist_route(c, dense) <= NEAR_NM:
            shore.append(c)
    log(f"{len(shore)} shoreline points")
    idx = [{"k": f"S{i + 1}", "kind": "shore", "lat": round(la, 5), "lon": round(lo, 5), "route_nm": round(dist_route((la, lo), dense), 2)}
           for i, (la, lo) in enumerate(sorted(shore, key=lambda p: -p[0]))]
    idx += [{"k": f"F{f['n']}", "kind": f["type"], "name": f.get("name", ""), "lat": f["lat"], "lon": f["lon"], "route_nm": f["route_nm"]}
            for f in feats if f["type"] not in ("DEPARE", "SEAARE", "LNDRGN", "CTNARE", "LAKARE", "RIVERS")]
    sd = OUT / "sheets"
    sd.mkdir(exist_ok=True)
    for old in sd.glob("*.jpg"):
        old.unlink()
    for s in range(0, len(idx), 25):
        sheet = Image.new("RGB", (5 * 244, 5 * 258), (255, 255, 255))
        d = ImageDraw.Draw(sheet)
        for k, o in enumerate(idx[s:s + 25]):
            x, y = (k % 5) * 244, (k // 5) * 258
            sheet.paste(chip(o["lat"], o["lon"]), (x + 2, y + 2))
            d.text((x + 4, y + 243), f"{o['k']} {o['kind'][:14]} {o['lat']:.4f},{o['lon']:.4f}", fill=(0, 0, 0))
        sheet.save(sd / f"sheet{s // 25 + 1:02d}.jpg", quality=78)
        o0, o1 = idx[s], idx[min(len(idx), s + 25) - 1]
        o0["sheet"] = f"sheet{s // 25 + 1:02d}.jpg"
        log(f"  sheet {s // 25 + 1}: {o0['k']} to {o1['k']}")
    for i, o in enumerate(idx):
        o["sheet"] = f"sheet{i // 25 + 1:02d}.jpg"
    (OUT / "index.json").write_text(json.dumps(idx, indent=0))
    (OUT / "features.json").write_text(json.dumps({
        "note": "Storm-shelter survey: chart (CHS ENC) and OpenStreetMap features within near_nm of the route between lat0 and lat1, "
                "with ~400 m Esri satellite chips on sheets/ (index.json). Built by scripts/shelter_survey.py. Not for navigation.",
        "lat0": LAT0, "lat1": LAT1, "near_nm": NEAR_NM, "generated": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()),
        "chart_object_counts": counts, "features": feats}, indent=0))
    log("done")


if __name__ == "__main__":
    main()
