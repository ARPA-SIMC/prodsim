#!/bin/bash

# Script for generating an orography database at the desired
# resolution on a desired grid

# define the area (e.g. Italy including minor islands)
XMIN=7.
XMAX=19.
YMIN=35.4
YMAX=47.1
NX=721 # 12*60+1
NY=703 # 11.7*60+1
GRIB_TEMPLATE=template.grib

# input points to a local directory on malina:/scratch/dcesari where
# globe dataset is stored and globalized for gdal
# vg6d_transform --trans-type=boxinter --sub-type=average --type=regular_ll \
vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll \
	       --nx=$NX --ny=$NY \
	       --x-min=$XMIN --y-min=$YMIN --x-max=$XMAX --y-max=$YMAX \
	       gdal,6.,35.,20.,48.:$SCRATCH/gis/raster/globe_30s/globe.vrt \
	       grib_api:${GRIB_TEMPLATE}:orog_tmp.grib

# replace missing values (~sea) with zeroes
vg6d_transform --trans-type=metamorphosis --sub-type=setinvalidto \
	       --maskbounds=0. orog_tmp.grib orog_hires.grib

rm -f orog_tmp.grib
