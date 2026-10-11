#!/usr/bin/env python3
"""Temporary: look for the YB race keys of the R2AK Proving Ground (Stage 1)."""
import json, urllib.request, urllib.error, sys
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
YEARS = [2016, 2017, 2018, 2019, 2022, 2023, 2024, 2025, 2026]
forms = []
for y in YEARS:
    yy = str(y)[2:]
    for f in ["r2akpg{y}", "r2akpg{yy}", "r2ak{y}pg", "r2ak{yy}pg", "pg{y}", "provingground{y}", "provingground{yy}", "r2akprovingground{y}",
              "r2ak{y}provingground", "r2akproving{y}", "r2ak{y}proving", "r2ak{y}leg1", "r2akleg1{y}", "r2akstage1{y}", "r2ak{yy}stage1",
              "r2ak{yy}s1", "r2ak{y}-pg", "r2ak{y}-stage1", "r2ak{y}_pg", "r2ak{y}prologue", "r2ak{y}a", "r2ak{y}1", "r2ak{y}s", "r2ak{y}p",
              "r2ak{y}stageone", "r2ak{y}ptv", "r2ak{y}pt", "r2akpt{y}", "r2ak{y}victoria", "r2ak{y}stage2", "r2ak{y}s2", "r2ak{yy}",
              "r2a{y}", "racetoalaska{y}", "rta{y}", "pgrace{y}", "r2ak{y}provinggrounds", "provinggrounds{y}"]:
        forms.append(f.format(y=y, yy=yy))
for k in dict.fromkeys(forms):
    try:
        raw = urllib.request.urlopen(urllib.request.Request(f"https://yb.tl/JSON/{k}/RaceSetup", headers=UA), timeout=20).read()
    except urllib.error.HTTPError as e:
        continue
    except Exception as e:
        print(k, "ERR", e, flush=True); continue
    try:
        s = json.JSONDecoder(strict=False).raw_decode(raw.decode("utf-8", "replace").strip())[0]
        print(k, "FOUND", repr(s.get("title")), len(s.get("teams") or []), s.get("start"), flush=True)
    except Exception:
        print(k, "nonjson", len(raw), raw[:80], flush=True)
