#!/usr/bin/env python3
"""Rebuild the Stage 2 south and inside legs of data/route.json from past human-powered tracks.

The route follows where past human-powered teams actually went: it is the cheapest water path
on the 220 m grid from build_route.py, where a cell costs less the more past tracks passed
within about 1 km of it. That keeps it on water (an average of coordinates would cut across
islands where teams chose different channels) while tracing the common line. It also passes
through Fitz Hugh Sound, inside Calvert Island.

Needs numpy, scipy, scikit-image, shapely and dl/land.npy written by build_route.py.
Run from the repository root: python3 scripts/build_consensus_route.py [path/to/land.npy]
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt
from skimage.draw import line
from skimage.graph import route_through_array
from shapely.geometry import LineString

ROOT = Path(__file__).resolve().parent.parent
X0, X1, Y0, Y1 = -131.9, -122.6, 48.0, 55.5
DY = 0.002; DX = DY / math.cos(math.radians(51.5))
land = np.load(sys.argv[1] if len(sys.argv) > 1 else "dl/land.npy")
H, W = land.shape


def rc(lat, lon): return (int(round((Y1 - lat) / DY)), int(round((lon - X0) / DX)))
def ll(r, c): return (Y1 - r * DY, X0 + c * DX)


# Lama Passage (past Bella Bella) is closed at 220 m: carve a corridor along the line past teams took
for (la0, lo0), (la1, lo1) in [((52.080, -128.040), (52.092, -128.110)), ((52.092, -128.110), (52.110, -128.137)), ((52.110, -128.137), (52.165, -128.140))]:
    for t in np.linspace(0, 1, 400):
        r, c = rc(la0 + (la1 - la0) * t, lo0 + (lo1 - lo0) * t)
        land[r - 1:r + 2, c - 1:c + 2] = False

# one boolean layer per track, widened to ~1 km, summed into a share of tracks per cell
tracks = []
for f in sorted((ROOT / "data" / "tracks").glob("r2ak*.json")):
    race = json.loads(f.read_text())
    for t in race["teams"]:
        pts = [(p[1], p[2]) for p in t["track"] if 48 < p[1] < 55.5]
        if len(pts) > 20:
            tracks.append((f"{race['key']} {t['name']}", pts))
count = np.zeros(land.shape, np.float32)
foot = np.ones((9, 9), bool)
for name, pts in tracks:
    m = np.zeros(land.shape, bool)
    for a, b in zip(pts, pts[1:]):
        (r0, c0), (r1, c1) = rc(*a), rc(*b)
        if max(abs(r1 - r0), abs(c1 - c0)) > 400:   # skip gaps over ~80 km (tracker outages)
            continue
        rr, cc = line(r0, c0, r1, c1)
        ok = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
        m[rr[ok], cc[ok]] = True
    count += binary_dilation(m, foot)
print(f"{len(tracks)} tracks", file=sys.stderr)
# share of the tracks that reach this far north, so the northern legs are not starved
reach = np.zeros(H, np.float32)
for name, pts in tracks:
    reach[: rc(max(p[0] for p in pts), 0)[0] + 1] += 1
share = count / np.maximum(reach, 1)[:, None]
dist = distance_transform_edt(~land) * DY * 60
cost = (1 + 3 * np.exp(-dist / 0.25)) / (0.08 + share)
cost[land] = 1e6


def route(vias):
    pts = []
    for a, b in zip(vias, vias[1:]):
        p, _ = route_through_array(cost, rc(*a), rc(*b), fully_connected=True, geometric=True)
        seg = [ll(r, c) for r, c in p]
        pts += seg if not pts else seg[1:]
    ls = LineString([(lo, la) for la, lo in pts])
    for tol in (0.004, 0.002, 0.001, 0.0005, 0.0):
        s = ls.simplify(tol, preserve_topology=False) if tol else ls
        if all(not land[rc(y0 + (y1 - y0) * t, x0 + (x1 - x0) * t)]
               for (x0, y0), (x1, y1) in zip(s.coords, list(s.coords)[1:])
               for t in np.linspace(0, 1, int(max(abs(x1 - x0) / DX, abs(y1 - y0) / DY)) + 2)):
            break
    return [[round(y, 4), round(x, 4)] for x, y in s.coords]


def nm(p): return sum(math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0]))) for a, b in zip(p, p[1:]))


VIC, SEYMOUR, FITZHUGH, BELLA, KTN = (48.4225, -123.372), (50.133, -125.35), (51.55, -127.86), (52.162, -128.14), (55.339, -131.637)
_, (ri, ci) = distance_transform_edt(land, return_indices=True)


def water(p):   # snap a via to the nearest water cell (narrows can be closed at 220 m)
    r, c = rc(*p)
    return ll(ri[r, c], ci[r, c])


VIC, SEYMOUR, FITZHUGH, BELLA, KTN = map(water, (VIC, SEYMOUR, FITZHUGH, BELLA, KTN))
south = route([VIC, SEYMOUR, BELLA])
inside = route([BELLA, KTN])
data = json.loads((ROOT / "data" / "route.json").read_text())
for leg in data["legs"]:
    if leg["key"] == "south":
        leg.update(points=south, name="Stage 2: Victoria to Bella Bella (where past human-powered teams went, inside Calvert Island)")
    if leg["key"] == "inside":
        leg.update(points=inside, name="Stage 2 north, inside: Bella Bella to Ketchikan (where past human-powered teams went)")
data["note"] = ("Stage 2 south and inside legs follow the common line of past human-powered tracks (build_consensus_route.py); "
                "Stage 1 and the outside leg are shortest-water tracks (build_route.py). For distances and ETAs, not navigation.")
(ROOT / "data" / "route.json").write_text(json.dumps(data, separators=(",", ":")))
print(f"south {nm(south):.1f} nm, {len(south)} pts; inside {nm(inside):.1f} nm, {len(inside)} pts", file=sys.stderr)
