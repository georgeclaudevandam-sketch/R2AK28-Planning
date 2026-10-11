#!/usr/bin/env python3
"""Tidal stream along the whole route from two bathymetry-driven ocean models, written to data/stream_field.json.

  WebTide, Northeast Pacific (ne_pac4)  DFO Bedford Institute / Institute of Ocean Sciences barotropic finite-element tidal
                                        model (Foreman et al.). Gives the u (east) and v (north) amplitude and Greenwich phase of
                                        eight constituents (M2 S2 N2 K2 K1 O1 P1 Q1) on a mesh fitted to the charted depths. Its
                                        mesh leaves out the Strait of Georgia and the Discovery Islands, and some narrow inlets.
  SalishSeaCast (UBC)                   3-D NEMO model of the Salish Sea run daily by UBC EOAS, on a 500 m grid. Hourly surface
                                        (0.5 m) velocities from its hindcast are fitted with UTide here (ten constituents and the
                                        mean, which carries the Fraser outflow and the estuarine surface flow). Covers Juan de Fuca,
                                        Haro Strait, the Strait of Georgia, Discovery Passage and the east end of Johnstone Strait.

Every mile along each route (Stage 1, and Stage 2 by the inside and the outside routes) the script takes SalishSeaCast where
it has water, else WebTide, and resolves the stream onto the route's heading there (+ = with you). It stores, per sample, the
along-track mean and the amplitude and phase of each constituent, plus for each constituent its nodal factor, frequency and
astronomical argument at the start of the 2028 window, so the planner can predict the stream at any time:

    along(t) = mean + sum_k f_k * A_k * cos(ph0_k + w_k * hours(t - t0) - g_k)

It also compares both models with the CHS/NOAA current stations (peak flood and ebb speeds over the window) and stores that.

Needs numpy, scipy, utide. Run from the repository root:
    python3 scripts/build_stream_field.py WEBTIDE_DATA_DIR      (the ne_pac4 folder of an installed WebTide)
Set SSC_CACHE to a folder to keep the SalishSeaCast downloads between runs.
"""
import csv, io, json, math, os, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
import utide
from utide._ut_constants import ut_constants
from utide.harmonics import FUV

ROOT = Path(__file__).resolve().parent.parent
WT = Path(sys.argv[1])
KN = 1.943844
STEP_NM = 1.0
CONST = ["M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1", "M4", "MS4"]
WT_CONST = ["M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1"]
ERDDAP = "https://salishsea.eos.ubc.ca/erddap/griddap"
SSC_U, SSC_V, SSC_GRID = "ubcSSg3DuGridFields1hV21-11", "ubcSSg3DvGridFields1hV21-11", "ubcSSnBathymetryV21-08"
FIT = ("2024-04-01T00:30:00Z", "2024-09-30T23:30:00Z")   # six months: separates K1/P1 and S2/K2
SSC_EDGE = 12            # skip grid cells this close to the model's open boundaries
CACHE = Path(os.environ.get("SSC_CACHE", ROOT / "dl" / "ssc"))


def log(*a): print(*a, file=sys.stderr, flush=True)


def nm(a, b):   # same haversine as index.html nmBetween
    R, r = 3440.065, math.pi / 180
    dl, dn = (b[0] - a[0]) * r, (b[1] - a[1]) * r
    x = math.sin(dl / 2) ** 2 + math.cos(a[0] * r) * math.cos(b[0] * r) * math.sin(dn / 2) ** 2
    return 2 * R * math.asin(math.sqrt(x))


class Chain:
    def __init__(s, pts):
        s.p = [tuple(p) for p in pts]
        s.cum = [0.0]
        for a, b in zip(s.p, s.p[1:]): s.cum.append(s.cum[-1] + nm(a, b))
        s.total = s.cum[-1]

    def at(s, a):
        a = max(0.0, min(s.total, a))
        i = max(1, next((k for k, c in enumerate(s.cum) if c >= a), len(s.cum) - 1))
        f = (a - s.cum[i - 1]) / ((s.cum[i] - s.cum[i - 1]) or 1)
        p, q = s.p[i - 1], s.p[i]
        return (p[0] + (q[0] - p[0]) * f, p[1] + (q[1] - p[1]) * f)

    def heading(s, a):   # radians from north, over a mile either side (as index.html trackHeading)
        p, q = s.at(a - 1), s.at(a + 1)
        return math.atan2((q[1] - p[1]) * math.cos(math.radians(p[0])), q[0] - p[0])


# ---------- WebTide
class WebTide:
    def __init__(s, d, box=(-138, -121, 46, 61)):
        n = np.loadtxt(d / "ne_pac4_ll.nod"); e = np.loadtxt(d / "ne_pac4.ele", dtype=int)[:, 1:] - 1
        s.lon, s.lat = n[:, 1], n[:, 2]
        ok = np.all((s.lon[e] > box[0]) & (s.lon[e] < box[1]) & (s.lat[e] > box[2]) & (s.lat[e] < box[3]), axis=1)
        s.e = e[ok]
        s.k = math.cos(math.radians(52))
        s.x, s.y = s.lon * s.k, s.lat
        s.tree = cKDTree(np.c_[s.x[s.e].mean(1), s.y[s.e].mean(1)])
        s.cu, s.cv = {}, {}
        for c in WT_CONST:
            a = np.loadtxt(d / f"{c}.barotropic.v2c", skiprows=3)
            idx = a[:, 0].astype(int) - 1
            cu = np.zeros(len(s.lon), complex); cv = np.zeros(len(s.lon), complex)
            cu[idx] = a[:, 1] * np.exp(-1j * np.radians(a[:, 2]))
            cv[idx] = a[:, 3] * np.exp(-1j * np.radians(a[:, 4]))
            s.cu[c], s.cv[c] = cu, cv

    def locate(s, lat, lon):
        px, py = lon * s.k, lat
        _, idx = s.tree.query([px, py], 12)
        for i in idx:
            a, b, c = s.e[i]
            x1, y1, x2, y2, x3, y3 = s.x[a], s.y[a], s.x[b], s.y[b], s.x[c], s.y[c]
            det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            if det == 0: continue
            l1 = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / det
            l2 = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / det
            l3 = 1 - l1 - l2
            if min(l1, l2, l3) >= -1e-9: return np.array([a, b, c]), np.array([l1, l2, l3])
        return None, None

    def uv(s, lat, lon):
        """{const: (Cu, Cv)} complex amplitudes in knots (A e^-ig), or None off the mesh"""
        nodes, w = s.locate(lat, lon)
        if nodes is None: return None
        return {"mean": (0j, 0j), **{c: (KN * (s.cu[c][nodes] * w).sum(), KN * (s.cv[c][nodes] * w).sum()) for c in WT_CONST}}


# ---------- SalishSeaCast
def erddap_csv(url, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=900) as r:
                rows = list(csv.reader(io.TextIOWrapper(r, "utf-8")))
            return rows[2:]
        except Exception as e:
            log(f"  retry {k + 1}: {e}")
            time.sleep(10 * (k + 1))
    raise RuntimeError(f"ERDDAP failed: {url}")


def q(s): return urllib.parse.quote(s, safe="(),:=")


class SalishSea:
    def __init__(s):
        CACHE.mkdir(parents=True, exist_ok=True)
        g = CACHE / "grid.npz"
        if not g.exists():
            rows = erddap_csv(f"{ERDDAP}/{SSC_GRID}.csv?" + q("longitude[0:1:897][0:1:397],latitude[0:1:897][0:1:397],bathymetry[0:1:897][0:1:397]"))
            a = np.array([[float(x) if x not in ("", "NaN") else np.nan for x in r] for r in rows])
            ny, nx = 898, 398
            np.savez(g, lon=a[:, 2].reshape(ny, nx), lat=a[:, 3].reshape(ny, nx), dep=np.nan_to_num(a[:, 4].reshape(ny, nx)))
        z = np.load(g)
        s.lon, s.lat, s.dep = z["lon"], z["lat"], z["dep"]
        ny, nx = s.lon.shape
        wet = s.dep > 0
        wet[:SSC_EDGE, :] = wet[-SSC_EDGE:, :] = False; wet[:, :SSC_EDGE] = wet[:, -SSC_EDGE:] = False
        yy, xx = np.nonzero(wet)
        s.cells = np.c_[yy, xx]
        s.tree = cKDTree(np.c_[s.lon[yy, xx] * math.cos(math.radians(49.5)), s.lat[yy, xx]])
        # angle of the grid's x axis from east, for rotating (u, v) on the grid to east/north
        dlon = np.gradient(s.lon, axis=1) * np.cos(np.radians(s.lat)); dlat = np.gradient(s.lat, axis=1)
        s.ang = np.arctan2(dlat, dlon)

    def cell(s, lat, lon, max_nm=0.4):
        d, i = s.tree.query([lon * math.cos(math.radians(49.5)), lat])
        if d * 60 > max_nm: return None
        return tuple(int(v) for v in s.cells[i])

    def series(s, y, x):
        f = CACHE / f"{y}_{x}.npz"
        if f.exists():
            z = np.load(f); return z["t"], z["u"], z["v"]
        t0, t1 = FIT
        ru = erddap_csv(f"{ERDDAP}/{SSC_U}.csv?" + q(f"uVelocity[({t0}):1:({t1})][(0.5000003)][({y}):1:({y})][({x - 1}):1:({x})]"))
        rv = erddap_csv(f"{ERDDAP}/{SSC_V}.csv?" + q(f"vVelocity[({t0}):1:({t1})][(0.5000003)][({y - 1}):1:({y})][({x}):1:({x})]"))
        def avg(rows):
            out = {}
            for r in rows: out.setdefault(r[0], []).append(float(r[4]) if r[4] not in ("", "NaN") else np.nan)
            ks = sorted(out)
            return np.array([datetime.fromisoformat(k.replace("Z", "+00:00")).timestamp() for k in ks]), np.array([np.nanmean(out[k]) for k in ks])
        tu, u = avg(ru); tv, v = avg(rv)
        n = min(len(tu), len(tv)); t, u, v = tu[:n], u[:n], v[:n]
        np.savez(f, t=t, u=u, v=v)
        return t, u, v

    def uv(s, lat, lon):
        c = s.cell(lat, lon)
        if c is None: return None
        y, x = c
        t, ug, vg = s.series(y, x)
        ok = np.isfinite(ug) & np.isfinite(vg)
        if ok.sum() < 24 * 120 or np.nanmax(np.abs(ug[ok])) == 0 and np.nanmax(np.abs(vg[ok])) == 0: return None
        a = s.ang[y, x]
        ue = (ug * math.cos(a) - vg * math.sin(a))[ok] * KN; vn = (ug * math.sin(a) + vg * math.cos(a))[ok] * KN
        tt = (t[ok] * 1e9).astype("int64").astype("datetime64[ns]")
        out = {}
        for name, comp in (("u", ue), ("v", vn)):
            co = utide.solve(tt, comp, lat=lat, method="ols", conf_int="none", trend=False, nodal=True, constit=CONST, verbose=False)
            out[name] = ({n: co["A"][i] * np.exp(-1j * np.radians(co["g"][i])) for i, n in enumerate(co["name"])}, co["mean"])
        res = {"mean": (complex(out["u"][1]), complex(out["v"][1]))}
        for c in CONST: res[c] = (out["u"][0].get(c, 0j), out["v"][0].get(c, 0j))
        return res


# ---------- astronomy for the 2028 window
def astro(t0, t1):
    names = list(ut_constants.const.name)
    lind = np.array([names.index(c) for c in CONST])
    dn = lambda d: d.timestamp() / 86400 + 719163.0     # datenum as utide uses it (1 = 1 Jan of year 1)
    tm = dn(t0 + (t1 - t0) / 2)
    F, U, V = FUV(np.array([dn(t0)]), dn(t0), lind, 50.0, [0, 0, 0, 0])
    Fm, _, _ = FUV(np.array([tm]), tm, lind, 50.0, [0, 0, 0, 0])
    freq = ut_constants.const.freq[lind]   # cycles per hour
    return {c: {"f": round(float(Fm[0, i]), 4), "ph0": round(float(((U[0, i] + V[0, i]) * 360) % 360), 3), "w": round(float(freq[i] * 360), 7)}
            for i, c in enumerate(CONST)}


def along(c, hd):
    """complex along-track amplitude from (Cu, Cv): component along heading hd (radians from north)"""
    return c[0] * math.sin(hd) + c[1] * math.cos(hd)


def predict(model, hd, t_hours, ast):
    x = np.full_like(t_hours, along(model["mean"], hd).real, dtype=float)
    for c in CONST:
        if c not in model: continue
        z = along(model[c], hd); A, g = abs(z), -math.degrees(math.atan2(z.imag, z.real))
        x += ast[c]["f"] * A * np.cos(np.radians(ast[c]["ph0"] + ast[c]["w"] * t_hours - g))
    return x


def main():
    route = json.loads((ROOT / "data" / "route.json").read_text())
    leg = {l["key"]: l["points"] for l in route["legs"]}
    chains = {"stage1": Chain(leg["stage1"]), "inside": Chain(leg["south"] + leg["inside"][1:]), "outside": Chain(leg["south"] + leg["outside"][1:])}
    cfg = json.loads((ROOT / "data" / "stations.json").read_text())
    t0 = datetime.fromisoformat(cfg["window"]["start"].replace("Z", "+00:00")); t1 = datetime.fromisoformat(cfg["window"]["end"].replace("Z", "+00:00"))
    ast = astro(t0, t1)
    wt = WebTide(WT)
    log("WebTide loaded")
    ssc = SalishSea()
    log("SalishSeaCast grid loaded")

    # every sample point, then the SalishSeaCast cells they need (fetched in parallel)
    samples = {k: [(a, ch.at(a), ch.heading(a)) for a in np.arange(0, ch.total + 1e-9, STEP_NM)] for k, ch in chains.items()}
    cells = sorted({c for ss in samples.values() for _, p, _ in ss for c in [ssc.cell(*p)] if c})
    log(f"{len(cells)} SalishSeaCast cells to fetch")
    def get(c):
        try: ssc.series(*c)
        except Exception as e: log(f"  cell {c} failed: {e}")
    with ThreadPoolExecutor(4) as ex:
        for i, _ in enumerate(ex.map(get, cells)):
            if i % 20 == 0: log(f"  {i}/{len(cells)}")

    memo = {}
    def model_at(p):
        c = ssc.cell(*p)
        if c is not None:
            if c not in memo:
                try: memo[c] = ssc.uv(*p)
                except Exception as e: log(f"  fit {c} failed: {e}"); memo[c] = None
            if memo[c] is not None: return 2, memo[c]
        m = wt.uv(*p)
        return (1, m) if m else (0, None)

    out = {}
    for k, ss in samples.items():
        rows = []
        for a, p, hd in ss:
            src, m = model_at(p)
            if not m: rows.append([0]); continue
            row = [src, round(along(m["mean"], hd).real, 3)]
            for c in CONST:
                z = along(m[c], hd) if c in m else 0j
                row += [round(abs(z), 3), round(-math.degrees(math.atan2(z.imag, z.real)) % 360)]
            rows.append(row)
        out[k] = rows
        n = [r[0] for r in rows]
        log(f"{k}: {len(rows)} samples, SalishSeaCast {n.count(2)}, WebTide {n.count(1)}, none {n.count(0)}")

    # check against the station tables: peak flood and ebb along each station's flood direction
    pred = json.loads((ROOT / "data" / "predictions.json").read_text())["stations"]
    hrs = np.arange(0, (t1 - t0).total_seconds() / 3600, 1 / 6)
    check = []
    for st in cfg["stations"]:
        sd = pred.get(st["key"])
        if st["kind"] != "current" or not sd or sd.get("flood_dir") is None: continue
        ev = [e for e in sd["events"] if e["type"] != "slack"]
        fl = np.percentile([e["kn"] for e in ev if e["type"] == "flood"], 90); eb = np.percentile([e["kn"] for e in ev if e["type"] == "ebb"], 90)
        row = {"key": st["key"], "table_flood": round(fl, 2), "table_ebb": round(eb, 2)}
        for name, src in (("ssc", 2), ("webtide", 1)):
            m = (ssc.uv(st["lat"], st["lon"]) if src == 2 and ssc.cell(st["lat"], st["lon"]) else None) if src == 2 else wt.uv(st["lat"], st["lon"])
            if not m: continue
            x = predict(m, math.radians(sd["flood_dir"]), hrs, ast)
            row[name + "_flood"] = round(float(np.percentile(x, 99.5)), 2); row[name + "_ebb"] = round(float(-np.percentile(x, 0.5)), 2)
        check.append(row); log(row)

    (ROOT / "data" / "stream_field.json").write_text(json.dumps({
        "note": "Tidal stream along the route from SalishSeaCast (UBC, src 2) and WebTide NE Pacific (DFO, src 1), resolved onto the route heading, + = with you, knots. "
                "Each sample (every step_nm along the chain): [src, mean, A1, g1, A2, g2, ...] in the order of const; [0] = no model. "
                "along(t) = mean + sum f*A*cos(ph0 + w*hours since t0 - g). Built by scripts/build_stream_field.py. For planning, not navigation.",
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"), "t0": cfg["window"]["start"], "step_nm": STEP_NM,
        "const": CONST, "astro": ast, "fit": {"salishsea": f"{SSC_U}/{SSC_V} surface 0.5 m, {FIT[0]} to {FIT[1]}", "webtide": "ne_pac4 barotropic"},
        "check": check, "chains": out}, separators=(",", ":")))
    log("wrote data/stream_field.json")


if __name__ == "__main__":
    main()
