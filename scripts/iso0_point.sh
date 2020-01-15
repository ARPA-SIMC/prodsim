#!/bin/bash

# shell script for retrieving height of 0 deg isotherm on selected points
#
# This shell script shows how to retrieve height of 0 deg isotherm
# from a gridded NWP field, interpolate it on some selected points
# with libsim tools, and obtain the result in a text file, e.g. for
# calculations related to radar products.

# get date from arguments
DATE="${1}-${2}-${3} ${4}:00"
DATEFILE=$1$2$3$4
# arkimet dataset
DS=http://arkimet.metarpa:8090/dataset/cosmo_5M_ita
# accumulation step
STEP='0 01'
LON_LIST=12.,10.
LAT_LIST=44.,46.

# define useful arkimet keys
# product: height of iso 0 (COSMO specific)
PROD_ISO0="GRIB1,80,201,84"
LEV_ISO0="GRIB1,4"

# clean old files
rm -f iso0.grib iso0_miss.grib

set -x
# split the query to avoid undesired fields because of "or" operator
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_ISO0; level: $LEV_ISO0" $DS > iso0.grib

vg6d_transform --trans-type=metamorphosis --sub-type=settoinvalid \
	       --maskbounds=-1000.,-998. iso0.grib iso0_miss.grib

vg6d_getpoint --lon="$LON_LIST" --lat="$LAT_LIST" \
	      --trans-type=inter --sub-type=near \
	      --output-format=grib_api_csv \
	      --output-keys=gacsv:simpleverdate,gacsv:lon,gacsv:lat,gacsv:value \
	      iso0_miss.grib iso0.csv

