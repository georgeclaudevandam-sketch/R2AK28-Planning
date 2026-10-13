#!/usr/bin/env python3
"""Take the zigzags out of every leg of data/route.json.

The consensus route follows the middle of past human-powered tracks every 0.2 nm, so where the tracks spread out it wanders
from side to side (a dogleg off Qualicum, a jog off Cape Lazo). Nobody rows that: they hold a heading and change it at a point.
So each leg is redrawn as straight runs. From each point the run goes as far along the route as a straight line can while
  - it stays on water at least a cell off the shore (data/water_mask.png, the 220 m shoreline grid; where that grid closes a narrows the old route went
    through, being within 150 m of the old route counts as water), and
  - it never strays more than MAX_OFF nm from the old route, nor leaves any of the old route more than MAX_OFF nm from it
    (so it keeps to the channel and the side of each island past teams used).
Run from the repository root after build_consensus_route.py: python3 scripts/smooth_route.py
"""
import json, math, sys
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
MAX_OFF = 0.5     # nm
NARROW = 0.08     # nm: a dry sample this close to the old route is a narrows the 220 m grid closes
STEP = 0.03       # nm between samples
from scipy.ndimage import binary_erosion
# water at least one 220 m cell from land, so a straight run keeps off the beach
MASK = binary_erosion(np.array(Image.open(ROOT / "data" / "water_mask.png").convert("L")) > 127)
Y1, DY, X0, DX = 55.5, 0.002, -131.9, 0.002 / math.cos(math.radians(51.5))


def wet(la, lo):
    j, i = int(round((Y1 - la) / DY)), int(round((lo - X0) / DX))
    return not (0 <= j < MASK.shape[0] and 0 <= i < MASK.shape[1]) or bool(MASK[j, i])


def nm(p): return sum(math.hypot((b[0] - a[0]) * 60, (b[1] - a[1]) * 60 * math.cos(math.radians(a[0]))) for a, b in zip(p, p[1:]))


def smooth(pts):
    k = 60 * math.cos(math.radians(sum(p[0] for p in pts) / len(pts)))
    xy = lambda p: (p[0] * 60, p[1] * k)
    # the old route, densely, with each sample's position along it
    dense, at = [], []
    for i, (a, b) in enumerate(zip(pts, pts[1:])):
        n = max(1, int(math.hypot(*np.subtract(xy(b), xy(a))) / STEP))
        for s in range(n):
            dense.append((a[0] + (b[0] - a[0]) * s / n, a[1] + (b[1] - a[1]) * s / n)); at.append(i)
    dense.append(tuple(pts[-1])); at.append(len(pts) - 1)
    D = np.array([xy(p) for p in dense])
    tree = cKDTree(D)

    def ok(i, j):
        a, b = np.array(xy(pts[i])), np.array(xy(pts[j]))
        L = np.linalg.norm(b - a); n = max(2, int(L / STEP))
        S = a + (b - a) * np.linspace(0, 1, n + 1)[:, None]
        d, _ = tree.query(S)
        if d.max() > MAX_OFF: return False
        for (x, y), dd in zip(S, d):
            if dd > NARROW and not wet(x / 60, y / k): return False
        # and the old route between i and j stays close to the new line
        seg = D[(np.array(at) >= i) & (np.array(at) < j)]
        if len(seg):
            t = np.clip(((seg - a) @ (b - a)) / max(L * L, 1e-12), 0, 1)
            if np.linalg.norm(seg - (a + t[:, None] * (b - a)), axis=1).max() > MAX_OFF: return False
        return True

    out, i = [pts[0]], 0
    while i < len(pts) - 1:
        j = i + 1
        # step out quickly, then back off to the furthest point that still works
        lo, hi = i + 1, i + 1
        while hi < len(pts) - 1 and ok(i, min(len(pts) - 1, hi * 2 - i)): lo = hi = min(len(pts) - 1, hi * 2 - i)
        hi = min(len(pts) - 1, hi * 2 - i)
        while hi - lo > 1:
            m = (lo + hi) // 2
            if ok(i, m): lo = m
            else: hi = m
        j = hi if ok(i, hi) else lo
        out.append(pts[j]); i = j
    return out


data = json.loads((ROOT / "data" / "route.json").read_text())
for leg in data["legs"]:
    old = leg["points"]
    new = smooth(old)
    leg["points"] = [[round(p[0], 5), round(p[1], 5)] for p in new]
    print(f"{leg['key']}: {len(old)} -> {len(new)} points, {nm(old):.1f} -> {nm(new):.1f} nm", file=sys.stderr)
data["note"] = data["note"].rstrip(".") + "; then straightened into runs that stay on water within 0.5 nm of that line (smooth_route.py)."
(ROOT / "data" / "route.json").write_text(json.dumps(data, separators=(",", ":")))
