#!/usr/bin/env python3
"""Import the R2AK 2015-2019 Swiftsure tracker archive into data/tracks.

The archive (from George's PC) is one GPX per year and stage, one <trk> per tracker, with a table of teams and results
(r2ak-tracks-2015-2019.csv, tracker_id = the GPX track name). Writes, in the same shape fetch_tracks.py writes for YB:
  - data/tracks/r2ak{year}-swiftsure.json: Stage 2 tracks of the human-powered teams (no sail), which the planner draws and
    in which it looks for the places teams stopped;
  - data/tracks/stage1/r2ak{year}-swiftsure.json: Stage 1 tracks of every boat.
The GPX <time> values are Pacific Daylight Time though marked Z, so 7 h is added. A finisher's track ends at its official
finish time; a boat that retired or was withdrawn ends at the first fix marked DNF or Withdrawn (the archive changes the status
once a day, so a finisher's "Finished" mark can come hours before the finish and is not used). 2015 teams already in r2ak2015-spot.json are left out of the 2015 file.

Run from the repository root: python3 scripts/import_swiftsure.py path/to/r2ak-baseline
"""
import csv, json, re, sys, xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_tracks import thin  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SRC = Path(sys.argv[1])
PDT = 7 * 3600
NS = {"g": "http://www.topografix.com/GPX/1/0"}
OUT = {"DNF", "Withdrawn"}


def epoch(s): return int(datetime.fromisoformat(s).timestamp()) if s else None


def boat_class(boat):
    b = boat.lower()
    if re.search(r"\bsup\b|paddleboard|airboard", b): return "SUP"
    if re.search(r"kayak|trak|passat|pygmy", b): return "Kayak"
    if "canoe" in b: return "Canoe"
    if "unlimited" in b: return "Human-powered"   # boat type not given
    return "Rowboat"


def human(row): return row["category"].startswith("Human") and not re.search("sail|pedal", row["boat"].lower())


def norm(n): return re.sub(r"[^a-z]", "", re.sub(r"^team\s+", "", n.lower()))


def tracks(path):
    out = {}
    for trk in ET.parse(path).getroot().findall("g:trk", NS):
        pts = []
        for p in trk.iter("{%s}trkpt" % NS["g"]):
            if p.findtext("g:desc", "", NS) in OUT: break
            t = int(datetime.strptime(p.findtext("g:time", "", NS), "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp()) + PDT
            pts.append((t, float(p.get("lat")), float(p.get("lon"))))
        out[trk.findtext("g:name", "", NS)] = sorted(pts)
    return out


rows = list(csv.DictReader(open(SRC / "r2ak-tracks-2015-2019.csv")))
spot = {norm(t["name"]) for t in json.loads((ROOT / "data/tracks/r2ak2015-spot.json").read_text())["teams"]}
for year in sorted({r["year"] for r in rows}):
    key = f"r2ak{year}-swiftsure"
    for stage in ("1", "2"):
        trk = tracks(SRC / "tracks" / f"r2ak-{year}-stage{stage}.gpx")
        rs = [r for r in rows if r["year"] == year and r["stage"] == stage]
        start = min(epoch(r["start_utc"]) for r in rs)
        teams = []
        for r in rs:
            fin = epoch(r["finish_utc"])
            pts = [p for p in trk.get(r["tracker_id"]) or [] if not fin or p[0] <= fin + 900]
            if stage == "2" and not human(r): continue
            if stage == "2" and year == "2015" and norm(r["team"]) in spot: continue
            pts = thin(pts, min_gap_s=300, min_nm=0.1) if stage == "2" else thin(pts)
            if len(pts) < 2: continue
            name = re.sub(r"^Team\s+", "", r["team"])
            if stage == "1":
                teams.append({"id": r["tracker_id"], "name": name, "boat": r["boat"] or None, "class": boat_class(r["boat"]) if human(r) else r["category"], "human": human(r), "status": r["status"], "elapsed_h": round((fin - epoch(r["start_utc"])) / 3600, 2) if fin else None, "track": pts})
            else:
                teams.append({"id": r["tracker_id"], "name": name, "boat": r["boat"] or None, "captain": r["skipper"] or None, "class": boat_class(r["boat"]),
                              "status": r["status"], "start": epoch(r["start_utc"]), "finished": fin, "elapsed_h": round((fin - epoch(r["start_utc"])) / 3600, 2) if fin else None, "track": pts})
        if stage == "1":
            out = {"key": key, "title": f"R2AK {year} Stage 1 (Swiftsure tracker archive)", "stage2_start": min(epoch(r["start_utc"]) for r in rows if r["year"] == year and r["stage"] == "2"), "teams": teams}
            (ROOT / "data/tracks/stage1" / f"{key}.json").write_text(json.dumps(out, separators=(",", ":")))
        else:
            out = {"key": key, "title": f"R2AK {year} (Swiftsure tracker archive)", "source": "Swiftsure tracker archive, R2AK 2015-2019 (copied from George's files, Oct 2026); times converted from PDT",
                   "note": "Human-powered teams only (rowboats, kayaks, canoes, paddleboards; no sail)." + (" Teams already in the 2015 SPOT file are left out." if year == "2015" else ""),
                   "start": start, "stop": max((p[0] for t in teams for p in t["track"]), default=start), "teams": teams}
            (ROOT / "data/tracks" / f"{key}.json").write_text(json.dumps(out, separators=(",", ":")))
        print(f"{key} stage {stage}: {len(teams)} teams, {sum(len(t['track']) for t in teams)} points" + ("" if stage == "1" else ": " + ", ".join(f"{t['name']} ({t['class']})" for t in teams)), file=sys.stderr)
