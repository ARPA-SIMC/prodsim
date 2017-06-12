#!/bin/bash

function lama_process() { # arguments reftime hour zdataset

echo "Processing $2"

# levels in LAMA dataset 
#LAMALEV=11,13,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40
#LAMALAY1=1,3,5,7,9,11,13,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40

LEVCONT=`seq 17 40`
ARKILAMALAY=`(for l in $LEVCONT; do echo -n "GRIB1,110,$l "; if [ $l -lt 40 ]; then echo -n "or ";fi; done)`
ARKILAMALEV=`(for l in $LEVCONT; do echo -n "GRIB1,109,$l "; if [ $l -lt 40 ]; then echo -n "or ";fi; done)`

ARKIDS="http://arkimet.metarpa:8090/dataset/lamaz"

# all variables at full levels (layers)
arki-query --data -o lama_fl.grib \
	   "Reftime:=$2; level: $ARKILAMALAY" \
	   $ARKIDS

# select only w al half levels
arki-query --data -o lama_hl.grib \
	   "Reftime:=$2; level: $ARKILAMALEV; product:GRIB1,,2,40;" \
	   $ARKIDS

# select T, U, V, P, Td Totprec at ~surface
arki-query --data -o lama_sl.grib \
	   "Reftime:=$2; level: GRIB1,105 or GRIB1,1; product:GRIB1,,2,11 or GRIB1,,2,33 or GRIB1,,2,34 or GRIB1,,2,1 or GRIB1,,2,17 or GRIB1,,2,61;" \
	   $ARKIDS

# filter the same layers from static dataset
arki-query --data -o lama_fl_const.grib \
	   "level: $ARKILAMALAY;" \
	   grib:$3

vg6d_transform --trans-type=vertint --sub-type=linear \
	       --trans-level-type=105,,105,105 \
	       lama_hl.grib lama_hl_to_fl.grib

cat lama_hl_to_fl.grib >> lama_fl.grib
rm -f lama_hl.grib lama_hl_to_fl.grib
mv lama_fl.grib lama_ua_$1.grib
mv lama_sl.grib lama_surf_$1.grib
mv lama_fl_const.grib lama_const_$1.grib

# do we need other upper air variables? B13003=rh, B11005=omega, B12102=td
#vg6d_transform --output-variable-list=B13003,B11005,B12102 \
#    lama_fl.grib lama_fl_newvar.grib

# convert all to first order, needed?
#/usr/libexec/ma_utils/grib_s2f.exe lama_ingv_test lama_ingv_f

}


ZDS=~eminguzzi/util/grib/lm_ope/LAMAZ_layers_20120606.grb
# levels in /home/eminguzzi/util/grib/lm_ope/LAMAZ_levels_20120606.grb

dt=${1::10}
while [ "$dt" -le "${2::10}" ]; do
    d=${dt::8}
    t=${dt:8:2}
    dtarki=`date -u --date "$d $t"  "+%Y-%m-%d %H:%M"`
    lama_process "$dt" "$dtarki" $ZDS
    dt=`date -u --date "$d $t 1 hour" "+%Y%m%d%H"`
done

