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
DS=http://arkimet.metarpa:8090/dataset/lmsmr6x54
# accumulation step
STEP='0 01'

# define useful arkimet keys
# product: wind, t and dew-point t, total precipitation, direct and diffuse radiation (COSMO specific)
PROD_CONST="GRIB1,,2,6 or GRIB1,,2,81"
PROD_WIND="GRIB1,,2,33 or GRIB1,,2,34"
PROD_TTD="GRIB1,,2,11 or GRIB1,,2,17"
PROD_PREC="GRIB1,,2,61"
PROD_RAD="GRIB1,,201,22 or GRIB1,,201,23"
# timerange: instantaneous, accumulated, averaged
TR_IST="GRIB1,0"
TR_ACC="GRIB1,4"
TR_AVG="GRIB1,3"
# level: height over surface (unspecified value), surface
LEV_HOS="GRIB1,105"
LEV_SURF="GRIB1,1"

# clean old files
rm -f uv.grib ttd.grib prec.grib rad.grib sd.grib trh.grib precacc.grib radavg.grib
# split the query to avoid undesired fields because of "or" operator
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_CONST; timerange: GRIB1,0,0; level: $LEV_SURF" $DS >> const.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_WIND; timerange: $TR_IST; level: $LEV_HOS" $DS >> uv.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_TTD; timerange: $TR_IST; level: $LEV_HOS" $DS >> ttd.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_PREC; timerange: $TR_ACC; level: $LEV_SURF" $DS  >> prec.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_RAD; timerange: $TR_AVG; level: $LEV_SURF" $DS  >> rad.grib

# compute height of surface (and keep land fraction)
vg6d_transform --output-variable-list=B10007,B29192 const.grib constz.grib
# compute wind speed and wind direction
vg6d_transform --output-variable-list=B11001,B11002 uv.grib sd.grib
# compute relative humidity (and keep temperature)
vg6d_transform --output-variable-list=B13003,B12101 ttd.grib trh.grib

# accumulate on desired interval
vg6d_transform --comp-stat-proc=1:1 --comp-step="$STEP" prec.grib precacc.grib
# average on desired interval
vg6d_transform --comp-stat-proc=0:0 --comp-step="$STEP" rad.grib radavg.grib
# ...
