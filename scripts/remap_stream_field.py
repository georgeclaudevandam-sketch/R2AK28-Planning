#!/usr/bin/env python3
"""Move data/stream_field.json onto a changed route without the models.

build_stream_field.py resolves the SalishSeaCast and WebTide streams onto the route's heading every mile, so a new route line
needs a new field. When the models can't be fetched (DFO took the WebTide installer down in October 2026), this rebuilds it
from what is already in the repository: data/stream_grid.json has the east/north stream (mean and the six main constituents)
on a half-mile grid within 3 nm of the old routes. Each new sample takes the grid cells within GRID_NM of it, weighted by
inverse distance, resolved onto the new heading; the four small constituents the grid lacks (K2, P1, Q1, MS4) and the model
source come from the nearest old sample. A sample stays with no model ([0]) where the nearest old sample had none.

Run from the repository root: python3 scripts/remap_stream_field.py path/to/old_route.json
"""
import base64, json, math, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_stream_field import Chain, STEP_NM, along  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
GRID_NM = 0.6
K = math.cos(math.radians(52))


def chains(route):
    leg = {l["key"]: l["points"] for l in route["legs"]}
    return {"stage1": Chain(leg["stage1"]), "inside": Chain(leg["south"] + leg["inside"][1:]), "outside": Chain(leg["south"] + leg["outside"][1:])}


def xy(p): return (p[1] * K * 60, p[0] * 60)


field = json.loads((ROOT / "data" / "stream_field.json").read_text())
G = json.loads((ROOT / "data" / "stream_grid.json").read_text())
old = chains(json.loads(Path(sys.argv[1]).read_text()))
new = chains(json.loads((ROOT / "data" / "route.json").read_text()))

g = G["grid"]
cells = np.frombuffer(base64.b64decode(G["cells"]), "<i4")
vals = np.frombuffer(base64.b64decode(G["vals"]), "<i2").reshape(len(cells), -1) / 100
glat = g["lat0"] + (cells // g["nLon"] + 0.5) * g["dlat"]; glon = g["lon0"] + (cells % g["nLon"] + 0.5) * g["dlon"]
gtree = cKDTree(np.c_[glon * K * 60, glat * 60])
gi = {c: i for i, c in enumerate(G["const"])}

out, moved = {}, []
for key, ch in new.items():
    o = old[key]
    orows = field["chains"][key]
    opts = [o.at(a) for a in np.arange(0, o.total + 1e-9, field["step_nm"])][: len(orows)]
    otree = cKDTree(np.array([xy(p) for p in opts]))
    rows = []
    for a in np.arange(0, ch.total + 1e-9, STEP_NM):
        p, hd = ch.at(a), ch.heading(a)
        _, j = otree.query(xy(p))
        base = orows[j]
        if not base[0]: rows.append([0]); continue
        idx = gtree.query_ball_point(xy(p), GRID_NM)
        if not idx: rows.append(base); continue
        d = np.hypot(*(np.array([xy((glat[i], glon[i])) for i in idx]) - xy(p)).T)
        w = 1 / np.maximum(d, 0.1) ** 2; w /= w.sum()
        v = (vals[idx] * w[:, None]).sum(0)
        row = [base[0], round(along((v[0], v[1]), hd), 3)]
        for n, c in enumerate(field["const"]):
            if c in gi:
                q = 2 + 4 * gi[c]
                z = along((complex(v[q], v[q + 1]), complex(v[q + 2], v[q + 3])), hd)
                row += [round(abs(z), 3), round(-math.degrees(math.atan2(z.imag, z.real)) % 360)]
            else:
                row += base[2 + 2 * n: 4 + 2 * n]
        rows.append(row)
    out[key] = rows
    print(f"{key}: {len(orows)} -> {len(rows)} samples, {sum(1 for r in rows if not r[0])} with no model", file=sys.stderr)

field["chains"] = out
field["note"] += " Remapped onto the current route from the stream grid by scripts/remap_stream_field.py."
(ROOT / "data" / "stream_field.json").write_text(json.dumps(field, separators=(",", ":")))
