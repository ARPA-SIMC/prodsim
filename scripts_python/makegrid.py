#!/usr/bin/python3

import subprocess

def convtoshp(name):
    try:
        subprocess.call(["ogr2ogr", name.replace("geojson", "shp"), name])
    except:
        print("conversion to shapefile failed")


intro = '''
{
  "type": "FeatureCollection",
  "features": [
'''

closing = '''
  ]
}
'''

rect = '''
    {
      "type": "Feature",
      "properties": {},
      "geometry": {
        "type": "Polygon",
        "coordinates": [
          [
            [
              %f,
              %f
            ],
            [
              %f,
              %f
            ],
            [
              %f,
              %f
            ],
            [
              %f,
              %f
            ],
            [
              %f,
              %f
            ]
          ]
        ]
      }
    }
'''

point = '''
    {
      "type": "Feature",
      "properties": {},
      "geometry": {
        "type": "Point",
        "coordinates": [
          %f,
          %f
        ]
      }
    }
'''



x0 = 6.
x1 = 19.
nx = 13
y0 = 35.
y1 = 49.
ny = 14

steps=(1,2,5,10)

for st in steps:
    jsname = f"box_{st}.geojson"
    js = open(jsname, "w")
    first = True
    js.write(intro)
    delta = 1./st
    for nj in range(ny*st):
        y = y0 + nj*delta
        for ni in range(nx*st):
            x = x0 + ni*delta
            coord = (x,y,x+delta,y,x+delta,y+delta,x,y+delta,x,y)
            if first: first = False
            else: js.write(",\n")
            js.write(rect % coord)
    js.write(closing)
    js.close()
    convtoshp(jsname)

    jsname = f"point_{st}.geojson"
    js = open(jsname, "w")
    first = True
    js.write(intro)
    delta = 1./st
    deltah = delta/2.
    for nj in range(ny*st):
        y = y0 + nj*delta + deltah
        for ni in range(nx*st):
            x = x0 + ni*delta + deltah
            coord = (x,y)
            if first: first = False
            else: js.write(",\n")
            js.write(point % coord)
    js.write(closing)
    js.close()
    convtoshp(jsname)

