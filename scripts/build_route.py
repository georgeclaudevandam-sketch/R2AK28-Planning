#!/usr/bin/env python3
"""Build data/route.json: shortest-water tracks for the R2AK course.

One-off helper, kept for reproducibility. Needs numpy, scipy, scikit-image, pillow,
shapely, plus two coastline files in dl/: gshhs_f.geojson (GSHHG full-resolution land
polygons for the route bbox, exported with basemap + basemap-data-hires) and
coastline.geojson (openwatersio/tide-database metadata/coastline.geojson, an
OpenStreetMap land-polygon extract).
"""
import json, math, numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt
from skimage.graph import route_through_array
from shapely.geometry import shape, LineString
X0,X1,Y0,Y1=-131.9,-122.6,48.0,55.5
DY=0.002; DX=DY/math.cos(math.radians(51.5))   # ~220 m cells
W=int((X1-X0)/DX); H=int((Y1-Y0)/DY)
img=Image.new('1',(W,H),0); d=ImageDraw.Draw(img)
def px(x,y): return ((x-X0)/DX,(Y1-y)/DY)
def drawgeom(g):
    for p in (g.geoms if hasattr(g,'geoms') else [g]):
        if p.geom_type!='Polygon': continue
        d.polygon([px(*c) for c in p.exterior.coords],fill=1)
        for i in p.interiors: d.polygon([px(*c) for c in i.coords],fill=0)
for fn in ('dl/gshhs_f.geojson','dl/coast.geojson'):
    for f in json.load(open(fn))['features']: drawgeom(shape(f['geometry']))
land=np.array(img,dtype=bool)
# Grenville Channel narrow portion is closed at GSHHS resolution: carve a corridor along its axis
for (la0,lo0),(la1,lo1) in [((53.395,-129.29),(53.605,-129.695)),((53.605,-129.695),(53.70,-129.85)),((53.70,-129.85),(53.80,-129.985)),((53.80,-129.985),(53.90,-130.09)),((53.90,-130.09),(53.96,-130.13)),((55.285,-131.545),(55.316,-131.603)),((55.316,-131.603),(55.336,-131.640))]:
    n=2000
    for t in np.linspace(0,1,n):
        r=int(round((Y1-(la0+(la1-la0)*t))/DY)); c=int(round(((lo0+(lo1-lo0)*t)-X0)/DX))
        land[r-1:r+2,c-1:c+2]=False
print('grid',W,H,land.mean())
dist=distance_transform_edt(~land)*DY*60   # nm-ish to nearest land
cost=np.where(land,np.inf,1+3*np.exp(-dist/0.25))
cost[land]=1e6
np.save('dl/land.npy',land)
def rc(lat,lon): return (int(round((Y1-lat)/DY)),int(round((lon-X0)/DX)))
def ll(r,c): return (Y1-r*DY, X0+c*DX)
def leg(a,b):
    p,_=route_through_array(cost,rc(*a),rc(*b),fully_connected=True,geometric=True)
    return [ll(r,c) for r,c in p]
def route(vias):
    pts=[]
    for a,b in zip(vias,vias[1:]):
        seg=leg(a,b); pts+= seg if not pts else seg[1:]
    ls=LineString([(lo,la) for la,lo in pts])
    # simplify but verify on raster
    for tol in (0.004,0.002,0.001,0.0005):
        s=ls.simplify(tol,preserve_topology=False)
        ok=True
        for (x0,y0),(x1,y1) in zip(s.coords,list(s.coords)[1:]):
            n=int(max(abs(x1-x0)/DX,abs(y1-y0)/DY))+2
            for t in np.linspace(0,1,n):
                r,c=rc(y0+(y1-y0)*t,x0+(x1-x0)*t)
                if land[r,c]: ok=False;break
            if not ok: break
        if ok: break
    print('tol',tol,len(s.coords))
    return [(round(y,4),round(x,4)) for x,y in s.coords]
def nm(pts):
    return sum(math.hypot((b[0]-a[0])*60,(b[1]-a[1])*60*math.cos(math.radians(a[0]))) for a,b in zip(pts,pts[1:]))
s1=route([(48.1175,-122.7480),(48.1550,-122.7400),(48.4040,-123.3900)])
south=route([(48.4225,-123.3720),(48.4050,-123.3900),(48.3870,-123.2950),(48.6000,-123.2050),(48.7550,-123.0700),(50.1330,-125.3500),(52.1620,-128.1400)])
outside=route([(52.1620,-128.1400),(55.0000,-131.1500),(55.2800,-131.5400),(55.3390,-131.6370)])
inside=route([(52.1620,-128.1400),(52.6000,-128.4600),(52.8700,-128.4900),(53.1000,-128.5400),(53.3300,-129.1500),(53.6000,-129.6830),(54.2000,-130.4000),(55.0000,-131.1500),(55.2800,-131.5400),(55.3390,-131.6370)])
for k,v in (('s1',s1),('south',south),('outside',outside),('inside',inside)): print(k,'%.1f nm'%nm(v))
json.dump({'s1':s1,'south':south,'outside':outside,'inside':inside},open('dl/autoroute.json','w'))
