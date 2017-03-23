#!/bin/bash

# define domain, orogsource and climsource
XMIN=7.
XMAX=19.
YMIN=35.4
YMAX=47.1
NX=721 # 12*60+1
NY=703 # 11.7*60+1
GRIB_TEMPLATE=template.grib
OROGSOURCE=/autofs/scratch-mod/dcesari/topo/globe_30s/globe.vrt
CLIMSOURCE=/autofs/scratch-mod/dcesari/climatology
SCRIPTS=/autofs/scratch-mod/ddalessandro/prodsim/scripts
SRC=/autofs/scratch-mod/ddalessandro/prodsim/src
INTERP=/autofs/scratch-mod/ddalessandro/prodsim/interpolamare

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
PROD_T="GRIB1,,2,11"
PROD_TD="GRIB1,,2,17"
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
rm -f const.grib uv.grib t.grib td.grib prec.grib rad.grib constz.grib sd.grib rh.grib precacc.grib radavg.grib
# split the query to avoid undesired fields because of "or" operator
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_CONST; timerange: GRIB1,0,0; level: $LEV_SURF" $DS > const.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_WIND; timerange: $TR_IST; level: $LEV_HOS" $DS > uv.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_T; timerange: $TR_IST; level: $LEV_HOS" $DS > t.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_TD; timerange: $TR_IST; level: $LEV_HOS" $DS > td.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_PREC; timerange: $TR_ACC; level: $LEV_SURF" $DS  > prec.grib
arki-query --data \
	   "reftime: ==$DATE; product: $PROD_RAD; timerange: $TR_AVG; level: $LEV_SURF" $DS > rad.grib

# compute height of surface (and throw away land fraction B29192)
vg6d_transform --output-variable-list=B10009 const.grib constz.grib
# compute wind speed and wind direction
vg6d_transform --output-variable-list=B11001,B11002 uv.grib sd.grib
# compute relative humidity
cat t.grib td.grib > ttd.grib
vg6d_transform --output-variable-list=B13003 ttd.grib rh.grib

# accumulate on desired interval
vg6d_transform --comp-stat-proc=1:1 --comp-step="$STEP" prec.grib precacc.grib
# average on desired interval
vg6d_transform --comp-stat-proc=0:0 --comp-step="$STEP" rad.grib radavg.grib
# ...


###########################  OROGPREPROC

# first make an identical transformation to grib
vg6d_transform --trans-type=none \
	       gdal,6.,35.,20.,48.:$OROGSOURCE \
	       grib_api:${GRIB_TEMPLATE}:orog_full.grib

# compute average, maximum and minimum for each cell
#for stat in average max min; do

    vg6d_transform --trans-type=boxinter --sub-type=average --type=regular_ll \
		   --nx=$NX --ny=$NY \
		   --x-min=$XMIN --y-min=$YMIN --x-max=$XMAX --y-max=$YMAX \
		   orog_full.grib orog_tmp.grib

# set to invalid sea points (only for gmrt, to avoid bathymetry)
    vg6d_transform --trans-type=metamorphosis --sub-type=settoinvalid \
		   --maskbounds=-15000.,0. orog_tmp.grib orog_tmp2.grib
# replace invalid values (~sea) with zeroes (for both)
    vg6d_transform --trans-type=metamorphosis --sub-type=setinvalidto \
		   --maskbounds=0. orog_tmp2.grib orog_hires_average.grib

    rm -f orog_tmp.grib orog_tmp2.grib
#done


###########################  CLIMPREPROC (to run just the first time)

# convert then cut in time and space
for file in $CLIMSOURCE/cru_v3_???_clim10.nc; do
    # local directory
    lfile=${file##*/}
    lfile=${lfile%.nc}
    # convert in grib setting correct parameter and unit
    # add setmisstonn to fill missing data
    cdo -f grb setparam,11.2 -addc,273.15 -setmisstonn $file $lfile.grib
    
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




##################################################

# climate orography interpolated over higher resolution grid
vg6d_transform --trans-type=inter --sub-type=bilin --output-format=grib_api:orog_hires_average.grib orog_cut.grib grib_api:orog_hires_average.grib:orog_average_cut_hires.grib 

# nwp orography interpolated over (the same) higher resolution grid
vg6d_transform --trans-type=inter --sub-type=bilin --output-format=grib_api:orog_hires_average.grib constz.grib grib_api:orog_hires_average.grib:constz_hires.grib 

# dati clima su grigliato finale
# vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll
#               --output-format=grib_api:orog_average_cut_hires.grib #
#                ${lfile}_1991_cut.grib
#               grib_api:orog_average_cut_hires.grib:${lfile}_interp.grib #


# clime interpolation over high resolution grid + temperature correction (moist adiabatic lapse rate) due to differences in orography 

for temp in tmn tmx tmp; do
vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib cru_v3_${temp}_clim10_1991_cut.grib grib_api:orog_average_cut_hires.grib:cru_v3_${temp}_clim10_interp.grib

$SRC/prodsim_vg6d_tcorr --tcorr-method=user --tgrad=-0.006 --input-orograhy=orog_average_cut_hires.grib --output-orograhy=orog_hires_average.grib cru_v3_${temp}_clim10_interp.grib cru_v3_${temp}_clim10_tcorr_saturo.grib

done 


# forecast interpolation over higher resolution grid

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib t.grib t_hires.grib 

#vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib sd.grib sd_hires.grib

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib uv.grib uv_hires.grib

#grib_filter $SCRIPTS/rules.txt $SCRIPTS/sd_hires.grib



#  temperature correction (moist adiabatic lapse rate)  due to differences in orography

$SRC/prodsim_vg6d_tcorr --tcorr-method=user --tgrad=-0.006 --input-orograhy=constz_hires.grib --output-orograhy=orog_hires_average.grib t_hires.grib previ_tcorr_saturo_oggi.grib

# compute maximum and minimum temperature fields

rm -f *_previ.grib
grib_filter $SCRIPTS/filtra_step.txt previ_tcorr_saturo_oggi.grib

b=24
c=48

for ((i=11, j=12; i<=23, j<=24; i+=1, j+=1)); do
    cdo max ${i}_previ.grib ${j}_previ.grib massime_oggi.grib
    if [ "$j" -ne "$b" ]; then
	grib_copy massime_oggi.grib ${j}_previ.grib
    fi
done

for ((i=35, j=36; i<=47, j<=48; i+=1, j+=1)); do
    cdo max ${i}_previ.grib ${j}_previ.grib massime_domani.grib
    if [ "$j" -ne "$c" ]; then
	grib_copy massime_domani.grib ${j}_previ.grib
    fi
done

rm -f *_previ.grib
grib_filter $SCRIPTS/filtra_step.txt previ_tcorr_saturo_oggi.grib

b=12
c=36

for ((i=0, j=1; i<=11, j<=12; i+=1, j+=1)); do
  
cdo min ${i}_previ.grib ${j}_previ.grib minime_oggi.grib

if [ "$j" -ne "$b" ]; then
grib_copy minime_oggi.grib ${j}_previ.grib
fi

done

for ((i=24, j=25; i<=35, j<=36; i+=1, j+=1)); do
  
cdo min ${i}_previ.grib ${j}_previ.grib minime_domani.grib

if [ "$j" -ne "$c" ]; then
grib_copy minime_domani.grib ${j}_previ.grib
fi

done

rm *_previ.grib

# compute tmax anomaly

$INTERP/anomalie_prova cru_v3_tmx_clim10_tcorr_saturo.grib massime_oggi.grib anomalie_oggi_tmx_clim10_corrette_Unipol.grib
$INTERP/anomalie_prova cru_v3_tmx_clim10_tcorr_saturo.grib massime_domani.grib anomalie_domani_tmx_clim10_corrette_Unipol.grib

cat anomalie_oggi_tmx_clim10_corrette_Unipol.grib anomalie_domani_tmx_clim10_corrette_Unipol.grib > anomalie_massime.grib

rm anomalie_oggi_tmx_clim10_corrette_Unipol.grib anomalie_domani_tmx_clim10_corrette_Unipol.grib

# compute tmin anomaly

$INTERP/anomalie_prova cru_v3_tmn_clim10_tcorr_saturo.grib minime_oggi.grib anomalie_oggi_tmn_clim10_corrette_Unipol.grib
$INTERP/anomalie_prova cru_v3_tmn_clim10_tcorr_saturo.grib minime_domani.grib anomalie_domani_tmn_clim10_corrette_Unipol.grib

cat anomalie_oggi_tmn_clim10_corrette_Unipol.grib anomalie_domani_tmn_clim10_corrette_Unipol.grib > anomalie_minime.grib

rm anomalie_oggi_tmn_clim10_corrette_Unipol.grib anomalie_domani_tmn_clim10_corrette_Unipol.grib




