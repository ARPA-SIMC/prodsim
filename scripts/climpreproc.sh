#!/bin/sh

XMIN=7.
XMAX=19.
YMIN=35.4
YMAX=47.1
NX=721 # 12*60+1
NY=703 # 11.7*60+1
GRIB_TEMPLATE=template.grib
OROGSOURCE=/autofs/scratch-mod/dcesari/topo/globe_30s/globe.vrt
CLIMSOURCE=/autofs/scratch-mod/dcesari/climatology

# Preprocessing for files like QESDI: Decadal surface meteorology
# climatologies from CRU TS3.0 data available at
# http://catalogue.ceda.ac.uk/uuid/432ac6151dcc4eff9f5874de3865784c

# convert then cut in time and space
for file in $CLIMSOURCE/cru_v3_???_clim10.nc; do
    # local directory
    lfile=${file##*/}
    lfile=${lfile%.nc}
    # convert in grib
    cdo -f grb copy $file $lfile.grib
    # keep only last decade 1991-2000
    grib_copy -w yearOfCentury=91 $lfile.grib ${lfile}_1991.grib
    vg6d_transform --trans-type=zoom --sub-type=coord \
		   --ilon=$XMIN --flon=$XMAX --ilat=$YMIN --flat=$YMAX \
		   ${lfile}_1991.grib ${lfile}_1991_cut.grib
done

# prepare corresponding orography (assuming it's compatible with the
# one used for the climatic dataset)

# first make an identical transformation to grib
vg6d_transform --trans-type=none \
	       gdal,6.,35.,20.,48.:$OROGSOURCE \
	       grib_api:${lfile}_1991_cut.grib:orog_full1.grib

# set to invalid sea points (only for gmrt, to avoid bathymetry)
vg6d_transform --trans-type=metamorphosis --sub-type=settoinvalid \
	       --maskbounds=-15000.,0. orog_full1.grib orog_full2.grib
# replace invalid values (~sea) with zeroes (for both)
vg6d_transform --trans-type=metamorphosis --sub-type=setinvalidto \
	       --maskbounds=0. orog_full2.grib orog_full3.grib

vg6d_transform --trans-type=boxinter --sub-type=average --type=regular_ll \
	       --output-format=grib_api:${lfile}_1991_cut.grib \
	       orog_full3.grib orog_cut.grib
