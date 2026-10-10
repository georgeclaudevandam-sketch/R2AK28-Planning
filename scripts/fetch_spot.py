#!/usr/bin/env python3
"""Fetch 2015 R2AK SPOT tracker logs from searunners.net and keep human-powered teams.

Source: http://searunners.net/data/R2AK_position_logs-sorted.tsv (name, obs_time
in Pacific Daylight Time, lat, lon) for every 2015 team. The file carries no boat
class, so data/tracks/spot2015-teams.json says which teams were human powered
(and their boat and result where known). Teams not listed there are logged so
the list can be completed.
"""
import json, math, sys, urllib.request
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
    raw = urllib.request.urlopen(urllib.request.Request(SRC, headers=UA), timeout=120).read().decode("utf-8", "replace")
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
    meta = json.loads(META.read_text()) if META.exists() else {}
    print(f"{len(teams)} teams in SPOT log", file=sys.stderr)
    out = []
    for name, pts in sorted(teams.items()):
        pts = sorted(set(pts))
        m = meta.get(name)
        print(f"  {'HP ' if m else '   '}{name!r:45} {len(pts):6} fixes  {datetime.fromtimestamp(pts[0][0], PDT):%m-%d %H:%M} -> "
              f"{datetime.fromtimestamp(pts[-1][0], PDT):%m-%d %H:%M}  last {pts[-1][1]:.2f},{pts[-1][2]:.2f}", file=sys.stderr)
        if not m:
            continue
        out.append({"id": name, "name": name.removeprefix("Team ").strip(), "boat": m.get("boat"), "captain": m.get("captain"),
                    "class": m.get("class", "Human powered"), "status": m.get("result"), "start": STAGE2_START,
                    "finished": None, "elapsed_h": m.get("elapsed_h"), "track": thin(pts)})
    missing = sorted(set(meta) - set(teams))
    if missing:
        print(f"listed but not in SPOT log: {missing}", file=sys.stderr)
    OUT.write_text(json.dumps({"key": "r2ak2015-spot", "title": "R2AK 2015 (SPOT trackers, searunners.net)",
                               "start": STAGE2_START, "source": SRC, "teams": out}, separators=(",", ":")))
    print(f"wrote {len(out)} human-powered teams", file=sys.stderr)


if __name__ == "__main__":
    main()
