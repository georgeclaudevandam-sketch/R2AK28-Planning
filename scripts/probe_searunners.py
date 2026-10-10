#!/usr/bin/env python3
"""Probe http://searunners.net/data/: crawl the directory listing and log what is there."""
import re, sys, urllib.request
from urllib.parse import urljoin
BASE = "http://searunners.net/data/"
UA = {"User-Agent": "r2ak28-planning (github.com/georgeclaudevandam-sketch/R2AK28-Planning)"}
seen, files = set(), []
def get(u, n=None):
    with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60) as r:
        return r.read(n) if n else r.read(), r.headers.get("Content-Type", "")
def crawl(u, depth=0):
    if u in seen or depth > 4 or len(seen) > 300: return
    seen.add(u)
    body, ct = get(u)
    links = re.findall(r'href="([^"?#]+)"', body.decode("utf-8", "replace"))
    for l in links:
        full = urljoin(u, l)
        if not full.startswith(BASE) or full == u or len(full) < len(u): continue
        if full.endswith("/"): crawl(full, depth + 1)
        else: files.append(full)
crawl(BASE)
print(f"{len(files)} files")
for f in files: print(f)
for f in files[:6] + [f for f in files if re.search(r'\.(csv|json|gpx|kml|txt|js)$', f, re.I)][:10]:
    try:
        b, ct = get(f, 1500); print("=====", f, ct); print(b.decode("utf-8", "replace")[:1500])
    except Exception as e: print("ERR", f, e)
