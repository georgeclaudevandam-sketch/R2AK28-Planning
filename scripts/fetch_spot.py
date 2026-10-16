#!/usr/bin/env python3
"""Fetch 2015 R2AK SPOT tracker logs from searunners.net and keep human-powered teams.

Source: http://searunners.net/data/R2AK_position_logs-sorted.tsv (name, obs_time
in Pacific Daylight Time, lat, lon) for every 2015 team. The file is the 28
teams that used human power at the start (oars, pedals or paddles; several also
carried sail). data/tracks/spot2015-teams.json gives each team's propulsion.
A team counts as finished at its first fix within 1.5 nm of Thomas Basin,
Ketchikan; other tracks are cut at their northernmost fix so the trip home
after retiring is left out.
"""
import json, math, sys, time, urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = "http://searunners.net/data/R2AK_position_logs-sorted.tsv"
META = ROOT / "data" / "tracks" / "spot2015-teams.json"
OUT = ROOT / "data" / "tracks" / "r2ak2015-spot.json"
PDT = timezone(timedelta(hours=-7))
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
STAGE2_START = int(datetime(2015, 6, 7, 12, 0, tzinfo=PDT).timestamp())


def thin(pts, min_gap_s=900, min_nm=0.3):
    out = []
    for p in pts:
        if out:
            q = out[-1]
            d = math.hypot((p[1] - q[1]) * 60, (p[2] - q[2]) * 60 * math.cos(math.radians(p[1])))
            if p[0] - q[0] < min_gap_s and d < min_nm:
                continue
        out.append(p)
    if pts and out[-1] != pts[-1]:
        out.append(pts[-1])
    return [[a, round(la, 4), round(lo, 4)] for a, la, lo in out]


def main():
    for attempt in range(4):
        try:
            raw = urllib.request.urlopen(urllib.request.Request(SRC, headers=UA), timeout=120).read().decode("utf-8", "replace")
            break
        except OSError as e:  # the archive host is slow at times
            if attempt == 3:
                raise
            print(f"retry {attempt + 1}: {e}", file=sys.stderr)
            time.sleep(15 * (attempt + 1))
    teams = {}
    for line in raw.splitlines()[1:]:
        f = [x.strip().strip('"') for x in line.split("\t")]
        if len(f) < 4:
            continue
        try:
            t = int(datetime.strptime(f[1], "%Y-%m-%d %H:%M:%S").replace(tzinfo=PDT).timestamp())
            la, lo = float(f[2]), float(f[3])
        except ValueError:
            continue
        if 40 < la < 60 and -140 < lo < -115:
            teams.setdefault(f[0], []).append((t, la, lo))
    meta = json.loads(META.read_text())["teams"] if META.exists() else {}
    print(f"{len(teams)} teams in SPOT log", file=sys.stderr)
    out = []
    for name, pts in sorted(teams.items()):
        pts = sorted(set(pts))
        m = meta.get(name)
        print(f"  {'HP ' if m else '   '}{name!r:45} {len(pts):6} fixes  {datetime.fromtimestamp(pts[0][0], PDT):%m-%d %H:%M} -> "
              f"{datetime.fromtimestamp(pts[-1][0], PDT):%m-%d %H:%M}  last {pts[-1][1]:.2f},{pts[-1][2]:.2f}", file=sys.stderr)
        if not m:
            continue
        pts = [p for p in pts if p[0] >= STAGE2_START - 3600]
        fin = next((p for p in pts if math.hypot((p[1] - 55.339) * 60, (p[2] + 131.637) * 60 * 0.57) < 1.5), None)
        if fin or m.get("official_elapsed_h"):
            if fin:
                pts = [p for p in pts if p[0] <= fin[0]]
            status = f"Finished {m['place']}" if m.get("place") else "Finished"
        else:
            north = max(range(len(pts)), key=lambda i: pts[i][1])
            pts = pts[:north + 1]
            status = "Did not finish"
        out.append({"id": name, "name": name.removeprefix("Team ").strip(), "boat": m.get("boat"), "captain": None,
                    "class": m.get("class"), "status": status, "start": STAGE2_START,
                    "finished": STAGE2_START + round(m["official_elapsed_h"] * 3600) if m.get("official_elapsed_h") else fin[0] if fin else None,
                    "elapsed_h": m.get("official_elapsed_h") or (round((fin[0] - STAGE2_START) / 3600, 2) if fin else None),
                    "place": m.get("place"), "track": thin(pts)})
    # Stage 1 (Port Townsend to Victoria, 4 June 2015): every team in the log, written like fetch_tracks.py's Stage 1 files
    s1 = []
    for name, pts in sorted(teams.items()):
        fx = sorted(p for p in set(pts) if STAGE2_START - 4 * 86400 < p[0] < STAGE2_START - 1800 and 47.9 < p[1] < 48.6 and -123.7 < p[2] < -122.5)
        if len(fx) < 3:
            continue
        m = meta.get(name)
        s1.append({"id": name, "name": name.removeprefix("Team ").strip(), "boat": (m or {}).get("boat"), "class": (m or {}).get("class") or "",
                   "human": bool(m), "track": thin(fx, min_gap_s=300, min_nm=0.15)})
    print(f"Stage 1: {len(s1)} teams with fixes, {sum(x['human'] for x in s1)} human-powered", file=sys.stderr)
    if s1:
        (OUT.parent / "stage1").mkdir(exist_ok=True)
        (OUT.parent / "stage1" / "r2ak2015-spot.json").write_text(json.dumps(
            {"key": "r2ak2015-spot", "title": "R2AK 2015 Stage 1 (SPOT)", "stage2_start": STAGE2_START, "teams": s1}, separators=(",", ":")))
    missing = sorted(set(meta) - set(teams))
    if missing:
        print(f"listed but not in SPOT log: {missing}", file=sys.stderr)
    OUT.write_text(json.dumps({"key": "r2ak2015-spot", "title": "R2AK 2015 (SPOT trackers via searunners.net)", "note": json.loads(META.read_text()).get("_note") + " Finish times are the official results (r2ak.com).",
                               "start": STAGE2_START, "source": SRC, "teams": out}, separators=(",", ":")))
    print(f"wrote {len(out)} teams", file=sys.stderr)


if __name__ == "__main__":
    main()
