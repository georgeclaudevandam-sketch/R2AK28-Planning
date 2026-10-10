#!/usr/bin/env python3
"""Water paths in to and out of each past campsite, written to data/camp_paths.json.

For each campsite and each Stage 2 route (inside and outside north of Bella Bella), the path in
leaves the route a little before the point abreast of the campsite and the path out rejoins it a
little after, so you row a loop rather than straight in and straight back. Both are least-cost
paths on the 220 m water grid from build_route.py that keep off the shore, the path out is pushed
away from the path in where there is room, and both are smoothed.

The campsites are the ones index.html finds in the past tracks (CAMPSITES), exported as a JSON list
of [lat, lon]. Paths are keyed by "lat,lon" to 4 decimals; a campsite with no path falls back to a
straight line in the planner.

Needs numpy, scipy, scikit-image, shapely and dl/land.npy.
Run from the repository root: python3 scripts/build_camp_paths.py campsites.json [path/to/land.npy]
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt
from skimage.graph import route_through_array
from shapely.geometry import LineString

ROOT = Path(__file__).resolve().parent.parent
X0, X1, Y0, Y1 = -131.9, -122.6, 48.0, 55.5
DY = 0.002; DX = DY / math.cos(math.radians(51.5))
camps = json.loads(Path(sys.argv[1]).read_text())
land = np.load(sys.argv[2] if len(sys.argv) > 2 else "dl/land.npy")
H, W = land.shape
dist = distance_transform_edt(~land) * DY * 60          # nm to the nearest land cell
_, (ri, ci) = distance_transform_edt(land, return_indices=True)


def rc(lat, lon): return (int(round((Y1 - lat) / DY)), int(round((lon - X0) / DX)))
def ll(r, c): return (Y1 - r * DY, X0 + c * DX)
def nm(a, b): return math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0])))
def length(p): return sum(nm(a, b) for a, b in zip(p, p[1:]))


def water(r, c): return (int(ri[r, c]), int(ci[r, c]))


def chain(points):
    cum = [0.0]
    for a, b in zip(points, points[1:]): cum.append(cum[-1] + nm(a, b))
    return points, cum


def at(ch, along):
    pts, cum = ch
    along = max(0, min(cum[-1], along))
    i = max(1, next((k for k, c in enumerate(cum) if c >= along), len(cum) - 1))
    f = (along - cum[i - 1]) / ((cum[i] - cum[i - 1]) or 1)
    a, b = pts[i - 1], pts[i]
    return (a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f)


def project(ch, p):
    pts, cum = ch
    i = min(range(len(pts)), key=lambda k: nm(pts[k], p))
    return cum[i], nm(pts[i], p)


def clear(p):   # every segment stays on water
    return all(not land[rc(a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)]
               for a, b in zip(p, p[1:]) for t in np.linspace(0, 1, int(max(abs(b[1] - a[1]) / DX, abs(b[0] - a[0]) / DY)) + 2))


def chaikin(p, n=3):
    for _ in range(n):
        q = [p[0]]
        for a, b in zip(p, p[1:]):
            q += [(0.75 * a[0] + 0.25 * b[0], 0.75 * a[1] + 0.25 * b[1]), (0.25 * a[0] + 0.75 * b[0], 0.25 * a[1] + 0.75 * b[1])]
        p = q + [p[-1]]
    return p


def tidy(cells):
    p = [ll(r, c) for r, c in cells]
    s = list(LineString([(b, a) for a, b in p]).simplify(0.0006).coords)
    s = [(y, x) for x, y in s]
    sm = chaikin(s)
    return sm if clear(sm) else s


def path(a, b, avoid=None):
    (r0, c0), (r1, c1) = a, b
    pad = 40
    R0, R1 = max(0, min(r0, r1) - pad), min(H, max(r0, r1) + pad + 1)
    C0, C1 = max(0, min(c0, c1) - pad), min(W, max(c0, c1) + pad + 1)
    cost = 1 + 4 * np.exp(-dist[R0:R1, C0:C1] / 0.12)
    if avoid is not None: cost = cost + 6 * avoid[R0:R1, C0:C1]
    cost[land[R0:R1, C0:C1]] = 1e6
    p, _ = route_through_array(cost, (r0 - R0, c0 - C0), (r1 - R0, c1 - C0), fully_connected=True, geometric=True)
    return [(r + R0, c + C0) for r, c in p]


route = json.loads((ROOT / "data" / "route.json").read_text())
leg = {l["key"]: [tuple(p) for p in l["points"]] for l in route["legs"]}
chains = {k: chain(leg["south"] + leg[k][1:]) for k in ("inside", "outside")}
out = {}
for site in camps:
    key = f"{site[0]:.4f},{site[1]:.4f}"
    cr = water(*rc(*site))
    for north, ch in chains.items():
        a, off = project(ch, site)
        if off > 12: continue
        d = max(0.6, min(3, 1.2 * off))          # leave and rejoin the route this far before and after
        A, B = at(ch, a - d), at(ch, a + d)
        pin = path(water(*rc(*A)), cr)
        near = np.zeros(land.shape, bool); rr, cc = zip(*pin); near[rr, cc] = True
        pout = path(cr, water(*rc(*B)), binary_dilation(near, np.ones((3, 3), bool)))
        pi, po = tidy(pin), tidy(pout)
        # a path far longer than the straight line means the grid closed a channel: leave that campsite to a straight line
        if max(length(pi), length(po)) > 2.5 * nm(site, A) + 1: continue
        rnd = lambda p: [[round(x, 4), round(y, 4)] for x, y in p]
        e = {"in": rnd(pi), "out": rnd(po), "nm_in": round(length(pi), 2), "nm_out": round(length(po), 2)}
        # south of Bella Bella the two routes are the same: store once
        if north == "outside" and key in out and out[key].get("inside") == e: continue
        out.setdefault(key, {})[north] = e
data = {"note": "Water paths in to and out of past campsites from the Stage 2 route (build_camp_paths.py). in starts on the route, out ends on it. Straight segments on a 220 m grid, smoothed; for planning, not navigation.",
        "sites": out}
(ROOT / "data" / "camp_paths.json").write_text(json.dumps(data, separators=(",", ":")))
print(f"{len(out)} campsites", file=sys.stderr)
