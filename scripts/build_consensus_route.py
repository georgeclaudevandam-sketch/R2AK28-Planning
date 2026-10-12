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

# Dodd Narrows (inside Gabriola Island, where 14 of the 15 past human-powered tracks went) is closed at 220 m as well: carve it
# along the line the 2026 teams' trackers recorded, from Stuart Channel to Nanaimo Harbour
DODD = [(49.1067, -123.8004), (49.1132, -123.8019), (49.1196, -123.8028), (49.1247, -123.8061), (49.1295, -123.8113), (49.1326, -123.8161),
        (49.1369, -123.8227), (49.1389, -123.8317), (49.1407, -123.8394), (49.1417, -123.8469), (49.1439, -123.8545)]
for (la0, lo0), (la1, lo1) in zip(DODD, DODD[1:]):
    for t in np.linspace(0, 1, 200):
        r, c = rc(la0 + (la1 - la0) * t, lo0 + (lo1 - lo0) * t)
        land[r - 1:r + 2, c - 1:c + 2] = False

# Enterprise Channel (inside Trial Islands, just out of Victoria), where every past human-powered track went, is closed at 220 m
# too: carve it along their common line
TRIAL = [(48.4030, -123.3290), (48.4045, -123.3185), (48.4070, -123.3100), (48.4082, -123.3020), (48.4095, -123.2955), (48.4140, -123.2915),
         (48.4185, -123.2885)]
for (la0, lo0), (la1, lo1) in zip(TRIAL, TRIAL[1:]):
    for t in np.linspace(0, 1, 200):
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
# Victoria Inner Harbour, in front of the Empress, out through the harbour entrance (too narrow for the 220 m grid; hand-drawn
# against the coastline). Stage 1 in build_route.py ends with the same path in reverse.
HARBOUR = [[48.4222, -123.3712], [48.4240, -123.3745], [48.4253, -123.3780], [48.4256, -123.3830], [48.4250, -123.3890], [48.4225, -123.3935], [48.4185, -123.3955]]
# then along the shore to Clover Point and through Enterprise Channel, on the common line of past human-powered tracks (a
# grid path wanders off it here: the shore is close and the channel is only open where carved)
OUT = [[48.4150, -123.3935], [48.4115, -123.3885], [48.4090, -123.3815], [48.4062, -123.3742], [48.4044, -123.3664], [48.4026, -123.3584],
       [48.4013, -123.3500], [48.4012, -123.3415], [48.4020, -123.3335]] + [list(p) for p in TRIAL]
k = next(k for k, q in enumerate(south) if q[1] > TRIAL[-1][1] + 0.002 and q[0] > TRIAL[-1][0])
south = HARBOUR + OUT + south[k:]
inside = route([BELLA, KTN])


# The grid path above keeps to water and to the common corridor, but where past tracks spread out (or hug a shore the 220 m grid
# closes) it can sit a mile or more off their middle. So pull it onto the tracks: every 0.2 nm along it, take each past track's
# nearest point (within 3 nm), keep the channel most of them used (the points within 0.5 nm of the medoid, the point closest to
# all the others), and move the route to their mean. Then smooth over about 1 nm.
from scipy.spatial import cKDTree
KX = 60 * math.cos(math.radians(51.5))
def xy(p): return (p[0] * 60, p[1] * KX)
trees = []
for name, pts in tracks:
    dense = []
    for a, b in zip(pts, pts[1:]):
        n = int(max(1, math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * KX) / 0.05))
        if n > 600:   # tracker gap
            continue
        dense += [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n)]
    if dense:
        trees.append((np.array(dense), cKDTree(np.array([xy(p) for p in dense]))))


def resample(p, step=0.2):
    out, carry = [p[0]], 0.0
    for a, b in zip(p, p[1:]):
        d = math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * KX)
        x = step - carry
        while x <= d:
            out.append((a[0] + (b[0] - a[0]) * x / d, a[1] + (b[1] - a[1]) * x / d)); x += step
        carry = d - (x - step)
    out.append(p[-1])
    return out


def wet(p): r, c = rc(*p); return not land[r, c]


def follow(p, keep=()):
    """keep: (lat, lon) boxes [la0, la1, lo0, lo1] left as they are (hand-drawn harbour)."""
    q = resample(p)
    out = []
    for s in q:
        if any(b[0] <= s[0] <= b[1] and b[2] <= s[1] <= b[3] for b in keep):
            out.append(s); continue
        near = []
        for dense, tree in trees:
            d, i = tree.query(xy(s), distance_upper_bound=3)
            if np.isfinite(d):
                near.append(dense[i])
        if len(near) < 3:
            out.append(s); continue
        P = np.array([xy(x) for x in near])
        D = np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :, 1])
        m = int(D.sum(1).argmin())
        grp = np.array(near)[D[m] < 0.5]
        mean = tuple(grp.mean(0))
        out.append(mean if wet(mean) or not wet(tuple(near[m])) else tuple(near[m]))
    # smooth over ~1 nm, but never onto land the unsmoothed point was off
    sm = []
    for i, s in enumerate(out):
        w = np.array(out[max(0, i - 2): i + 3]).mean(0)
        sm.append(tuple(w) if wet(tuple(w)) or not wet(s) else s)
    sm[0], sm[-1] = out[0], out[-1]
    ls = LineString([(lo, la) for la, lo in sm]).simplify(0.0003, preserve_topology=False)
    return [[round(y, 4), round(x, 4)] for x, y in ls.coords]


nh = len(HARBOUR)
south = HARBOUR + follow(south[nh - 1:], keep=[(48.40, 48.43, -123.40, -123.36)])[1:]
inside = follow(inside)
# South of Klemtu the 2026 teams split: most came up the west side (notes, Lillian Signed Up To Suffer, Laughing Roomba), two
# crossed from the east along 52.55 N. The medoid flips between the two groups and leaves a dogleg east of everyone's line, so
# hand-draw the west group's common line from 52.505 N to the mouth of Klemtu Passage.
KLEMTU = [[52.5097, -128.4742], [52.5208, -128.4803], [52.5311, -128.4847], [52.5398, -128.4898], [52.5479, -128.4947], [52.5539, -128.4975]]
a = next(i for i, p in enumerate(inside) if p[0] > 52.505)
b = next(i for i, p in enumerate(inside) if p[0] > 52.560)
inside = inside[:a] + KLEMTU + inside[b:]
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
