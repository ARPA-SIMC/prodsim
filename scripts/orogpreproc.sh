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


# OROGSOURCE points to a shared directory in ARPA-SIMC LAN
# from https://www.ngdc.noaa.gov/mgg/topo/globe.html
OROGSOURCE=/autofs/scratch-mod/dcesari/topo/globe_30s/globe.vrt
# from http://www.marine-geo.org/portals/gmrt/about.php
# OROGSOURCE=/autofs/scratch-mod/dcesari/topo/gmrt/GMRTv3_3_20170216topo.tif

# first make an identical transformation to grib
vg6d_transform --trans-type=none \
	       gdal,6.,35.,20.,48.:$OROGSOURCE \
	       grib_api:${GRIB_TEMPLATE}:orog_full.grib

# compute average, maximum and minimum for each cell
for stat in average max min; do

    vg6d_transform --trans-type=boxinter --sub-type=$stat --type=regular_ll \
		   --nx=$NX --ny=$NY \
		   --x-min=$XMIN --y-min=$YMIN --x-max=$XMAX --y-max=$YMAX \
		   orog_full.grib orog_tmp.grib

# set to invalid sea points (only for gmrt, to avoid bathymetry)
    vg6d_transform --trans-type=metamorphosis --sub-type=settoinvalid \
		   --maskbounds=-15000.,0. orog_tmp.grib orog_tmp2.grib
# replace invalid values (~sea) with zeroes (for both)
    vg6d_transform --trans-type=metamorphosis --sub-type=setinvalidto \
		   --maskbounds=0. orog_tmp2.grib orog_hires_$stat.grib

    rm -f orog_tmp.grib orog_tmp2.grib
done
