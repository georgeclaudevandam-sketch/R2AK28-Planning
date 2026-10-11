#!/usr/bin/env python3
"""Probe: can the CHS and NOAA chart services return LIGHTS features (characteristic, period, colour, sectors) for a box?"""
import json, sys, urllib.parse, urllib.request
UA = {"User-Agent": "r2ak28-planning"}
SVC = {"chs": "https://egisp.dfo-mpo.gc.ca/arcgis/rest/services/chs/ENC_MaritimeChartService/MapServer/exts/MaritimeChartService/MapServer",
       "noaa": "https://gis.charttools.noaa.gov/arcgis/rest/services/MCS/NOAAChartDisplay/MapServer/exts/MaritimeChartService/MapServer"}
BOX = {"chs": (-128.56, 52.62, -128.48, 52.67), "noaa": (-131.70, 55.30, -131.60, 55.36)}
def get(u):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=90) as r: return r.read().decode("utf-8", "replace")
    except Exception as e: return f"ERR {e}"
for k, svc in SVC.items():
    print("=====", k)
    print(get(svc + "?f=json")[:1500])
    x0, y0, x1, y1 = BOX[k]
    for layers in ("all", "all:0,1,2,3,4,5,6,7"):
        q = {"geometry": json.dumps({"xmin": x0, "ymin": y0, "xmax": x1, "ymax": y1, "spatialReference": {"wkid": 4326}}), "geometryType": "esriGeometryEnvelope",
             "sr": "4326", "layers": layers, "tolerance": "2", "mapExtent": f"{x0},{y0},{x1},{y1}", "imageDisplay": "800,600,96", "returnGeometry": "true", "f": "json"}
        r = get(svc + "/identify?" + urllib.parse.urlencode(q))
        print("--- identify", layers, len(r)); print(r[:4000])
    print("--- layers"); print(get(svc + "/layers?f=json")[:2000])
