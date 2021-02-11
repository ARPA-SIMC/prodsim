#!/bin/bash

# shell script for preliminary preprocessing of NWP data
#
# This shell script shows how to retrieve some typical NWP grib
# products from arkimet archive and preprocess them with libsim in
# order to obtain data more suitable for the preparation of end-user
# meteorological products.
#
# Final processed parameters are temperature, relative humidity, wind
# speed and direction, total precipitation, direct and diffuse
# radiation plus constant orography and fraction of land.

# get date from arguments
DATE="${1}-${2}-${3} ${4}:00"
DATEFILE=$1$2$3$4
# arkimet dataset
DS=http://arkimet.metarpa:8090/dataset/cosmo_5M_ita
# accumulation step
STEP='0 01'

# define useful arkimet keys
# product: wind, t and dew-point t, total precipitation, direct and diffuse radiation (COSMO specific)
PROD_ORO="GRIB1,,2,6"
PROD_FRLAND="GRIB1,,2,81"
PROD_WIND="GRIB1,,2,33 or GRIB1,,2,34"
PROD_T="GRIB1,,2,11"
PROD_TD="GRIB1,,2,17"
PROD_TS="GRIB1,,2,85"
PROD_PREC="GRIB1,,2,61"
PROD_SNOW="GRIB1,,2,78 or GRIB1,,2,79"
PROD_RAD="GRIB1,,201,22 or GRIB1,,201,23"
# timerange: instantaneous, accumulated, averaged
TR_IST="GRIB1,0"
TR_ACC="GRIB1,4"
TR_AVG="GRIB1,3"
# level: height over surface (unspecified value), surface
LEV_HOS="GRIB1,105"
LEV_SURF="GRIB1,1"
LEV_SOIL="GRIB1,111"

# clean old files
rm -f oro.grib frland.grib uv.grib t.grib td.grib ts.grib prec.grib rad.grib constz.grib sd.grib rh.grib precacc.grib radavg.grib

set -x
# split the query to avoid undesired fields because of "or" operator
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_ORO; timerange: GRIB1,0,0; level: $LEV_SURF" $DS > oro.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_FRLAND; timerange: GRIB1,0,0; level: $LEV_SURF" $DS > frland.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_WIND; timerange: $TR_IST; level: $LEV_HOS" $DS > uv.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_T; timerange: $TR_IST; level: $LEV_HOS" $DS > t.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_TD; timerange: $TR_IST; level: $LEV_HOS" $DS > td.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_TS; timerange: $TR_IST; level: $LEV_SOIL" $DS > ts.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_PREC or $PROD_SNOW; timerange: $TR_ACC; level: $LEV_SURF" $DS  > prec.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_RAD; timerange: $TR_AVG; level: $LEV_SURF" $DS > rad.grib

# compute height of surface
vg6d_transform --output-variable-list=B10009 oro.grib oroz.grib
# compute wind speed and wind direction
vg6d_transform --output-variable-list=B11001,B11002 uv.grib sd.grib
# compute relative humidity
cat t.grib td.grib > ttd.grib
vg6d_transform --output-variable-list=B13003 ttd.grib rh.grib

# accumulate on desired interval
vg6d_transform --comp-stat-proc=1:1 --comp-step="$STEP" prec.grib precacc.grib
# compute total snow (and keep total precipitation)
vg6d_transform --output-variable-list=B13011,B13201 precacc.grib precacctot.grib
# average on desired interval
vg6d_transform --comp-stat-proc=0:0 --comp-step="$STEP" rad.grib radavg.grib

# keep temperature and dew point temperature only on land points
vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.5,1.5 \
	       --coord-format=grib_api --coord-file=frland.grib \
	       ttd.grib ttd_landonly.grib

