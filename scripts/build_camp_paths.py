#!/usr/bin/env python3
"""Water paths in to and out of each past campsite, written to data/camp_paths.json.

For each campsite and each Stage 2 route (inside and outside north of Bella Bella), the path in
leaves the route before the point abreast of the campsite and the path out rejoins it after, so you
row a loop rather than straight in and straight back. Both follow least-cost paths on the 220 m
water grid from build_route.py that keep off the shore, the path out pushed away from the path in
where there is room. Where they meet the route they peel off and merge in smoothly: each starts
along the route's heading and curves (a cubic Bezier) onto the water path, and the place it leaves
the route and the length of the curve are chosen, among those that stay on water, for the gentlest
turns and the least extra rowing. No sharp turn off the route to a stop.

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
from skimage.graph import route_through_array, MCP_Geometric
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


def log(*a): print(*a, file=sys.stderr)


def hd(a, b): return math.degrees(math.atan2((b[1] - a[1]) * math.cos(math.radians(a[0])), b[0] - a[0]))


def point_along(p, s):
    acc = 0.0
    for a, b in zip(p, p[1:]):
        d = nm(a, b)
        if acc + d >= s: f = (s - acc) / (d or 1); return (a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f)
        acc += d
    return p[-1]


def after(p, s):   # the part of polyline p beyond distance s
    acc = 0.0
    for k, (a, b) in enumerate(zip(p, p[1:])):
        acc += nm(a, b)
        if acc > s: return p[k + 1:]
    return [p[-1]]


def max_turn(p, skip_end=0.15, step=0.1):
    """sharpest change of heading, in degrees per 0.1 nm, ignoring the last bit into the campsite"""
    L = length(p)
    q = [point_along(p, x) for x in np.arange(0, max(step, L - skip_end), step)]
    h = [hd(a, b) for a, b in zip(q, q[1:]) if nm(a, b) > 1e-4]
    return max([abs((b - a + 180) % 360 - 180) for a, b in zip(h, h[1:])] or [0])


def bezier(a, ta, b, tb, n=24):
    """cubic from a leaving along heading ta to b arriving along heading tb (degrees), handles 0.4 of the chord"""
    h = 0.4 * nm(a, b) / 60
    k = math.cos(math.radians(a[0]))
    c1 = (a[0] + h * math.cos(math.radians(ta)), a[1] + h * math.sin(math.radians(ta)) / k)
    c2 = (b[0] - h * math.cos(math.radians(tb)), b[1] - h * math.sin(math.radians(tb)) / k)
    return [((1 - t) ** 3 * a[0] + 3 * (1 - t) ** 2 * t * c1[0] + 3 * (1 - t) * t * t * c2[0] + t ** 3 * b[0],
             (1 - t) ** 3 * a[1] + 3 * (1 - t) ** 2 * t * c1[1] + 3 * (1 - t) * t * t * c2[1] + t ** 3 * b[1]) for t in np.linspace(0, 1, n)]


class Tree:
    """least-cost water paths from one cell (the campsite) to anywhere in a window around it"""
    def __init__(s, src, box, avoid=None):
        s.R0, s.R1, s.C0, s.C1 = box
        cost = 1 + 4 * np.exp(-dist[s.R0:s.R1, s.C0:s.C1] / 0.12)
        if avoid is not None: cost = cost + 6 * avoid[s.R0:s.R1, s.C0:s.C1]
        cost[land[s.R0:s.R1, s.C0:s.C1]] = 1e6   # not impassable: a campsite in a pond on the grid still gets a path
        s.m = MCP_Geometric(cost, fully_connected=True)
        s.m.find_costs([(src[0] - s.R0, src[1] - s.C0)])

    def to(s, cell):   # cells from cell to the source
        r, c = cell
        if not (s.R0 <= r < s.R1 and s.C0 <= c < s.C1): return None
        try: p = s.m.traceback((r - s.R0, c - s.C0))
        except Exception: return None
        return [(a + s.R0, b + s.C0) for a, b in p][::-1]


def rough(cells):
    p = [ll(r, c) for r, c in cells]
    s = list(LineString([(b, a) for a, b in p]).simplify(0.0006).coords)
    return [(y, x) for x, y in s]


def peel(ch, tree, site, a_join, sign):
    """water path between the route at a_join and the campsite, its route end tangent to the route.
    sign -1: leaving the route (path runs route -> campsite); +1: rejoining (built backwards, then reversed)"""
    A = at(ch, a_join)
    th = hd(at(ch, a_join - 0.15), at(ch, a_join + 0.15)) + (0 if sign < 0 else 180)   # heading away from the route end
    cells = tree.to(water(*rc(*A)))
    if not cells: return None
    g = rough(cells)                 # route end -> campsite (its nearest water)
    A = g[0]
    best, n = None, length(g)
    for L in sorted({0.25, 0.4, 0.6, 0.9, 1.3, 1.8} | {round(f * n, 2) for f in (0.3, 0.5, 0.7)}):
        if L < 0.05 or L > n - 0.05: continue
        P, P2 = point_along(g, L), point_along(g, L + 0.1)
        b = bezier(A, th, P, hd(P, P2))
        p = b + after(g, L)
        sm = chaikin(p[len(b) - 2:], 2)
        p = b[:-2] + sm if clear(b[:-2] + sm) else p
        if not clear(p): continue
        t, n = max_turn(p), length(p)
        score = n * (1 + max(0.0, t - 20) / 20)
        if best is None or score < best[0]: best = (score, p, t)
    if best is None:   # no curve stays on water: the water path as it is (smoothed where that stays on water)
        sm = chaikin(g); p = sm if clear(sm) else g
        best = (n * 10, p, max_turn(p))   # last resort: any curved candidate at another join point wins over this
    return best


route = json.loads((ROOT / "data" / "route.json").read_text())
leg = {l["key"]: [tuple(p) for p in l["points"]] for l in route["legs"]}
chains = {k: chain(leg["south"] + leg[k][1:]) for k in ("inside", "outside")}
out = {}
worst = []
for site in camps:
    key = f"{site[0]:.4f},{site[1]:.4f}"
    cr = water(*rc(*site))
    for north, ch in chains.items():
        a, off = project(ch, site)
        if off > 12: continue
        lo, hi = at(ch, max(0, a - 6)), at(ch, a + 6)
        rr = [rc(*p)[0] for p in (lo, hi, site)]; cc = [rc(*p)[1] for p in (lo, hi, site)]
        box = (max(0, min(rr) - 40), min(H, max(rr) + 41), max(0, min(cc) - 40), min(W, max(cc) + 41))
        tin = Tree(cr, box)
        # path in: leave the route somewhere from abreast of the campsite to 5 nm before it
        cand = [peel(ch, tin, site, a - d, -1) for d in np.arange(max(0.4, min(2.5, 0.8 * off)), max(0.4, min(2.5, 0.8 * off)) + 4.01, 0.25) if a - d >= 0]
        cand = [c for c in cand if c]
        if not cand: log(key, north, "no clear path in"); continue
        _, pi, ti = min(cand, key=lambda c: c[0])
        near = np.zeros(land.shape, bool)
        for x in np.linspace(0, length(pi), int(length(pi) / 0.05) + 2): near[rc(*point_along(pi, x))] = True
        tout = Tree(cr, box, binary_dilation(near, np.ones((5, 5), bool)))
        cand = [peel(ch, tout, site, a + d, +1) for d in np.arange(max(0.4, min(2.5, 0.8 * off)), max(0.4, min(2.5, 0.8 * off)) + 4.01, 0.25)]
        cand = [c for c in cand if c]
        if not cand: log(key, north, "no clear path out"); continue
        _, po, to = min(cand, key=lambda c: c[0])
        po = po[::-1]
        # a path far longer than the straight line means the grid closed a channel: leave that campsite to a straight line
        if length(pi) > 2.5 * nm(site, pi[0]) + 1.5 or length(po) > 2.5 * nm(site, po[-1]) + 1.5: log(key, north, f"too long {length(pi):.1f} {length(po):.1f} off {off:.1f}"); continue
        worst.append((max(ti, to), key, north))
        rnd = lambda p: [[round(x, 4), round(y, 4)] for x, y in p]
        e = {"in": rnd(pi), "out": rnd(po), "nm_in": round(length(pi), 2), "nm_out": round(length(po), 2)}
        # south of Bella Bella the two routes are the same: store once
        if north == "outside" and key in out and out[key].get("inside") == e: continue
        out.setdefault(key, {})[north] = e
worst.sort(reverse=True)
print("sharpest turns (deg per 0.1 nm):", [(round(t), k, n) for t, k, n in worst[:6]], file=sys.stderr)
data = {"note": "Water paths in to and out of past campsites from the Stage 2 route (build_camp_paths.py). in starts on the route, out ends on it. Straight segments on a 220 m grid, smoothed; for planning, not navigation.",
        "sites": out}
(ROOT / "data" / "camp_paths.json").write_text(json.dumps(data, separators=(",", ":")))
print(f"{len(out)} campsites", file=sys.stderr)
