#!/usr/bin/env python3
"""Fetch past R2AK tracks of human-powered teams from YB Tracking.

For every candidate YB race key that exists, reads RaceSetup (teams, boat
classes, finish times) and AllPositions3 (packed binary positions), keeps the
human-powered teams, thins their tracks and writes data/tracks/<key>.json plus
an index at data/tracks/index.json.

Human-powered = tagged Human Powered, Kayak, Rowboat, SUP, Canoe, Pedal or
Pedal Power by the race organisers in YB.
"""
from __future__ import annotations

import json
import math
import struct
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "tracks"
YB = "https://yb.tl"
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
HUMAN_TAGS = {"human powered", "kayak", "rowboat", "row boat", "sup", "canoe", "pedal", "pedal power", "paddle", "rowing"}
YEARS = [2015, 2016, 2017, 2018, 2019, 2022, 2023, 2024, 2025, 2026]
SUFFIXES = ["", "s1", "stage1"]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url: str, tries: int = 4, timeout: int = 30):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 410, 500):
                return None
            if e.code in (429, 503):
                time.sleep(20 * (i + 1))
            log(f"  HTTP {e.code} retry {i + 1}: {url}")
        except Exception as e:
            log(f"  {type(e).__name__} {e} retry {i + 1}: {url}")
        time.sleep(2 * (i + 1))
    return None


def decode_positions(buf: bytes) -> dict[int, list[tuple[int, float, float]]]:
    """Decode YB AllPositions3: big-endian, newest fix first, later fixes as deltas."""
    flags = buf[0]
    has_alt, has_dtf, has_laps, has_pc = flags & 1, flags & 2, flags & 4, flags & 8
    ref = struct.unpack_from(">I", buf, 1)[0]
    off, teams = 5, {}
    while off + 4 <= len(buf):
        tid, n = struct.unpack_from(">HH", buf, off)
        off += 4
        pts, prev = [], None
        for _ in range(n):
            if buf[off] & 0x80:
                w, dlat, dlon = struct.unpack_from(">Hhh", buf, off)
                off += 6
                if has_alt:
                    off += 2
                if has_dtf:
                    off += 2
                    if has_laps:
                        off += 1
                if has_pc:
                    off += 2
                at, lat, lon = prev[0] - (w & 0x7FFF), prev[1] + dlat, prev[2] + dlon
            else:
                t, lat, lon = struct.unpack_from(">Iii", buf, off)
                off += 12
                if has_alt:
                    off += 2
                if has_dtf:
                    off += 4
                    if has_laps:
                        off += 1
                if has_pc:
                    off += 4
                at = ref + t
            prev = (at, lat, lon)
            pts.append(prev)
        teams[tid] = sorted((a, la / 1e5, lo / 1e5) for a, la, lo in pts)
    return teams


def thin(pts, min_gap_s=900, min_nm=0.3):
    """Keep a fix when 15 min or 0.3 nm has passed since the last kept one."""
    out = []
    for p in pts:
        if not (40 < p[1] < 60 and -140 < p[2] < -115):
            continue
        if out:
            q = out[-1]
            d = math.hypot((p[1] - q[1]) * 60, (p[2] - q[2]) * 60 * math.cos(math.radians(p[1])))
            if p[0] - q[0] < min_gap_s and d < min_nm:
                continue
        out.append(p)
    if pts and out and out[-1] != pts[-1] and 40 < pts[-1][1] < 60:
        out.append(pts[-1])
    return [[int(a), round(la, 4), round(lo, 4)] for a, la, lo in out]


def race(key: str):
    log(f"{key}: RaceSetup")
    raw = get(f"{YB}/JSON/{key}/RaceSetup")
    if not raw:
        return None
    try:
        setup = json.JSONDecoder(strict=False).raw_decode(raw.decode("utf-8", "replace").strip())[0]
    except ValueError as e:
        log(f"  JSON error {e}")
        log(f"  not JSON ({len(raw)} bytes): {raw[:300]!r}")
        return None
    if not isinstance(setup, dict) or not setup.get("teams"):
        log(f"  no teams; keys {list(setup)[:20] if isinstance(setup, dict) else type(setup).__name__}")
        return None
    tags = {t["id"]: t["name"] for t in setup.get("tags", [])}
    human_ids = {i for i, n in tags.items() if n.strip().lower() in HUMAN_TAGS}
    log(f"{key}: '{setup.get('title')}' {len(setup['teams'])} teams, tags {sorted(tags.values())}")
    teams = [t for t in setup["teams"] if human_ids & set(t.get("tags") or [])]
    log(f"  {len(teams)} human-powered teams")
    if not teams:
        return {"key": key, "title": setup.get("title"), "start": setup.get("start"), "teams": []}
    buf = get(f"{YB}/BIN/{key}/AllPositions3", tries=3, timeout=120)
    if not buf:
        log("  no AllPositions3")
        return None
    pos = decode_positions(buf)
    allpts = [p for v in pos.values() for p in v]
    if allpts:
        log(f"  decoded {len(pos)} tracks, {len(allpts)} fixes, lat {min(p[1] for p in allpts):.2f}..{max(p[1] for p in allpts):.2f}, "
            f"lon {min(p[2] for p in allpts):.2f}..{max(p[2] for p in allpts):.2f}")
    out = []
    for t in teams:
        pts = thin(pos.get(t["id"], []))
        cls = [tags[i] for i in t.get("tags") or [] if i in tags and i in human_ids and tags[i].lower() != "human powered"]
        start = t.get("start") or setup.get("start")
        fin = t.get("finishedAt")
        out.append({
            "id": t["id"], "name": t.get("name"), "boat": t.get("model"), "captain": t.get("captain"),
            "class": cls[0] if cls else "Human powered", "status": t.get("explain") or t.get("status"),
            "start": start, "finished": fin, "elapsed_h": round((fin - start) / 3600, 2) if fin and start else None,
            "track": pts,
        })
        log(f"  {t.get('name')!r:40} {len(pts):5} pts  {'finished %.1f d' % ((fin - start) / 86400) if fin and start else t.get('explain') or t.get('status')}")
    return {"key": key, "title": setup.get("title"), "start": setup.get("start"), "stop": setup.get("stop"), "teams": out}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    keys = sys.argv[1:] or [f"r2ak{y}{s}" for s in SUFFIXES for y in reversed(YEARS)] + [f"r2ak{str(y)[2:]}" for y in YEARS]
    index = []
    for key in dict.fromkeys(keys):
        try:
            r = race(key)
        except Exception as e:  # keep going with other years
            log(f"{key}: FAILED {type(e).__name__}: {e}")
            continue
        if not r or not r["teams"]:
            continue
        (OUT / f"{key}.json").write_text(json.dumps(r, separators=(",", ":")))
        index.append({"key": key, "title": r["title"], "start": r["start"], "teams": len(r["teams"]),
                      "finished": sum(1 for t in r["teams"] if t["finished"])})
        (OUT / "index.json").write_text(json.dumps({"source": "YB Tracking (yb.tl)", "races": index}, indent=1))
    (OUT / "index.json").write_text(json.dumps({"source": "YB Tracking (yb.tl)", "races": index}, indent=1))
    log(f"wrote {len(index)} races")


if __name__ == "__main__":
    main()
