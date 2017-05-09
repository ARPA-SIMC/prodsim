#!/bin/bash

# define domain, orogsource and climsource
XMIN=6.
XMAX=18.5
YMIN=35.9
YMAX=47.1
NX=751 # 12.5*60+1
NY=673 # 11.2*60+1
GRIB_TEMPLATE=template.grib
OROGSOURCE=/autofs/scratch-mod/dcesari/topo/globe_30s/globe.vrt
CLIMSOURCE=/autofs/scratch-mod/dcesari/climatology
SCRIPTS=/autofs/scratch-mod/ddalessandro/prodsim/scripts
SRC=/autofs/scratch-mod/ddalessandro/prodsim/src
INTERP=/autofs/scratch-mod/ddalessandro/prodsim/interpolamare
DS=http://arkimet.metarpa:8090/dataset/lmsmr6x54

# accumulation step
STEP='0 01'

# define useful arkimet keys
# product: wind, t and dew-point t, total precipitation, direct and diffuse radiation (COSMO specific)
PROD_CONST="GRIB1,,2,6 or GRIB1,,2,81"
PROD_WIND="GRIB1,,2,33 or GRIB1,,2,34"
PROD_T="GRIB1,,2,11"
PROD_TD="GRIB1,,2,17"
#PROD_TS="GRIB1,,2,85"
PROD_PREC="GRIB1,,2,61"
PROD_RAD="GRIB1,,201,22 or GRIB1,,201,23"
PROD_SDI="GRIB1,,201,141"
PROD_MASK="GRIB1,,2,81"
# timerange: instantaneous, accumulated, averaged
TR_IST="GRIB1,0"
TR_ACC="GRIB1,4"
#TR_ACC="GRIB1,4,,12h or GRIB1,4,,24h or GRIB1,4,,36h or GRIB1,4,,48h"
TR_AVG="GRIB1,3"
# level: height over surface (unspecified value), surface
LEV_HOS="GRIB1,105"
LEV_SURF="GRIB1,1"
LEV_SOIL="GRIB1,111"

# clean old files
rm -f const.grib uv.grib t.grib td.grib prec.grib rad.grib constz.grib sd.grib rh.grib precacc.grib radavg.grib

date +"%H-%M"
MM=$(date +"%M") # minute
hh=$(date +"%H") # hour
gg=$(date +"%d") # day
mm=$(date +"%m") # month
yyyy=$(date +"%Y") # year

if [ "$hh" -lt 20 ]; then # for run00

arki-query --data \
	   "reftime: =today 00:00; product: $PROD_CONST; timerange: GRIB1,0,0; level: $LEV_SURF" $DS > const.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_WIND; timerange: $TR_IST; level: $LEV_HOS" $DS > uv.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_T; timerange: $TR_IST; level: $LEV_HOS" $DS > t.grib # COSMO I2 447x532
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_TD; timerange: $TR_IST; level: $LEV_HOS" $DS > td.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_TS; timerange: $TR_IST; level: $LEV_SOIL" $DS > ts.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_PREC; timerange: $TR_ACC; level: $LEV_SURF" $DS  > prec.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_RAD; timerange: $TR_AVG; level: $LEV_SURF" $DS > rad.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_SDI; timerange: $TR_IST" $DS > sdi.grib
arki-query --data \
	   "reftime: =today 00:00; product: $PROD_MASK" $DS > mask.grib
fi

if [ "$hh" -ge 19 ]; then # for run12


arki-query --data \
	   "reftime: =today 12:00; product: $PROD_CONST; timerange: GRIB1,0,0; level: $LEV_SURF" $DS > const.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_WIND; timerange: $TR_IST; level: $LEV_HOS" $DS > uv.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_T; timerange: $TR_IST; level: $LEV_HOS" $DS > t.grib # COSMO I2 447x532
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_TD; timerange: $TR_IST; level: $LEV_HOS" $DS > td.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_TS; timerange: $TR_IST; level: $LEV_SOIL" $DS > ts.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_PREC; timerange: $TR_ACC; level: $LEV_SURF" $DS  > prec.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_RAD; timerange: $TR_AVG; level: $LEV_SURF" $DS > rad.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_SDI; timerange: $TR_IST" $DS > sdi.grib
arki-query --data \
	   "reftime: =today 12:00; product: $PROD_MASK" $DS > mask.grib
fi


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
	       grib_api:$SCRIPTS/${GRIB_TEMPLATE}:orog_full.grib
                                       # orog_full.grib 1680x1560
# compute average, maximum and minimum for each cell
#for stat in average max min; do

    vg6d_transform --trans-type=boxinter --sub-type=average --type=regular_ll \
		   --nx=$NX --ny=$NY \
		   --x-min=$XMIN --y-min=$YMIN --x-max=$XMAX --y-max=$YMAX \
		   orog_full.grib orog_tmp.grib

# set to invalid sea points (only for gmrt, to avoid bathymetry)
    vg6d_transform --trans-type=metamorphosis --sub-type=settoinvalid \
		   --maskbounds=-15000.,0. orog_tmp.grib orog_tmp2.grib  #maskbounds defines the range of values to become invalid
# replace invalid values (~sea) with zeroes (for both)
    vg6d_transform --trans-type=metamorphosis --sub-type=setinvalidto \
		   --maskbounds=0. orog_tmp2.grib orog_hires_average.grib  #maskbounds sets the constant value to be used
                                                # orog_hires_average.grib 751x673
    rm -f orog_tmp.grib orog_tmp2.grib
#done




###########################  CLIMPREPROC (to run just the first time)

## convert then cut in time and space
#for file in $CLIMSOURCE/cru_v3_???_clim10.nc; do
#    # local directory
#    lfile=${file##*/}
#    lfile=${lfile%.nc}
#    # convert in grib setting correct parameter and unit
#    # add setmisstonn to fill missing data
#    cdo -f grb setparam,11.2 -addc,273.15 -setmisstonn $file $lfile.grib
   
#    # keep only last decade 1991-2000
#    grib_copy -w yearOfCentury=91 $lfile.grib ${lfile}_1991.grib
#    vg6d_transform --trans-type=zoom --sub-type=coord \
#		   --ilon=$XMIN --flon=$XMAX --ilat=$YMIN --flat=$YMAX \
#		   ${lfile}_1991.grib ${lfile}_1991_cut.grib
#done
# prepare corresponding orography (assuming it's compatible with the
# one used for the climatic dataset)

# first make an identical transformation to grib
#vg6d_transform --trans-type=none \
#	       gdal,6.,35.,20.,48.:$OROGSOURCE \
#	       grib_api:${lfile}_1991_cut.grib:orog_full1.grib
## set to invalid sea points (only for gmrt, to avoid bathymetry)
#vg6d_transform --trans-type=metamorphosis --sub-type=settoinvalid \
#	       --maskbounds=-15000.,0. orog_full1.grib orog_full2.grib
## replace invalid values (~sea) with zeroes (for both)
#vg6d_transform --trans-type=metamorphosis --sub-type=setinvalidto \
#	       --maskbounds=0. orog_full2.grib orog_full3.grib
#vg6d_transform --trans-type=boxinter --sub-type=average --type=regular_ll \
#	       --output-format=grib_api:${lfile}_1991_cut.grib \
#	       orog_full3.grib orog_cut.grib
#                             # orog_cut.grib 25x25 


##################################################

# climate orography interpolated over higher resolution grid
vg6d_transform --trans-type=inter --sub-type=bilin --output-format=grib_api:orog_hires_average.grib $SCRIPTS/orog_cut.grib grib_api:orog_hires_average.grib:orog_average_cut_hires.grib 
# orog_average_cut_hires.grib 751x673

# nwp orography interpolated over (the same) higher resolution grid
vg6d_transform --trans-type=inter --sub-type=bilin --output-format=grib_api:orog_hires_average.grib constz.grib grib_api:orog_hires_average.grib:constz_hires.grib 
# constz_hires.grib 751x673

# dati clima su grigliato finale
# vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll
#               --output-format=grib_api:orog_average_cut_hires.grib #
#                ${lfile}_1991_cut.grib
#               grib_api:orog_average_cut_hires.grib:${lfile}_interp.grib #


# clime interpolation over high resolution grid + temperature correction (moist adiabatic lapse rate) due to differences in orography 

for temp in tmn tmx tmp; do
vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib $SCRIPTS/cru_v3_${temp}_clim10_1991_cut.grib grib_api:orog_average_cut_hires.grib:cru_v3_${temp}_clim10_interp.grib

$SRC/prodsim_vg6d_tcorr --tcorr-method=user --tgrad=-0.006 --input-orograhy=orog_average_cut_hires.grib --output-orograhy=orog_hires_average.grib cru_v3_${temp}_clim10_interp.grib cru_v3_${temp}_clim10_tcorr_saturo.grib

done 

# mask interpolation over higher resolution grid

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll \
	       --output-format=grib_api:orog_average_cut_hires.grib \
	       mask.grib mask_hires.grib

# forecast interpolation over higher resolution grid

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib t.grib t_hires.grib 

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib sd.grib sd_hires.grib

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib uv.grib uv_hires.grib

vg6d_transform --trans-type=inter --sub-type=bilin --type=regular_ll --output-format=grib_api:orog_average_cut_hires.grib precacc.grib precacc_hires.grib

grib_filter $SCRIPTS/rules.txt $SCRIPTS/sd_hires.grib

grib_copy $SCRIPTS/32_sd.grib wind_speed.grib



#  temperature correction (moist adiabatic lapse rate)  due to differences in orography

$SRC/prodsim_vg6d_tcorr --tcorr-method=user --tgrad=-0.006 --input-orograhy=constz_hires.grib --output-orograhy=orog_hires_average.grib t_hires.grib previ_tcorr_saturo_oggi.grib

# compute maximum and minimum temperature fields

rm -f *_previ.grib
grib_filter $SCRIPTS/filtra_step.txt previ_tcorr_saturo_oggi.grib

rm -f anomalieTemperatureMassime_*.grib anomalieTemperatureMinime_*.grib precipitazioni_*.grib temperatura_*.grib vento_*.grib


# number of bands BB

BB3=$(grib_count precacc_hires.grib)
BB4=$(grib_count previ_tcorr_saturo_oggi.grib)
BB5=$(grib_count wind_speed.grib)

printf -v BB3 "%02d" $BB3
printf -v BB4 "%02d" $BB4
printf -v BB5 "%02d" $BB5


# NN

NN3=$(grib_get -w count=1 -p endStep precacc_hires.grib)

t5=$(grib_get -w count=1 -p stepRange previ_tcorr_saturo_oggi.grib)
t6=$(grib_get -w count=2 -p stepRange previ_tcorr_saturo_oggi.grib)
NN4=$((t6-t5))

t7=$(grib_get -w count=1 -p stepRange wind_speed.grib)
t8=$(grib_get -w count=2 -p stepRange wind_speed.grib)
NN5=$((t8-t7))

printf -v NN3 "%02d" $NN3
printf -v NN4 "%02d" $NN4
printf -v NN5 "%02d" $NN5

###########################################################################################à
#run00

if [ "$hh" -lt 20 ]; then

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

# number of bands BB?
BB1=$(grib_count anomalie_massime.grib)
BB2=$(grib_count anomalie_minime.grib)
printf -v BB1 "%02d" $BB1
printf -v BB2 "%02d" $BB2

t1=$(grib_get -w count=1 -p stepRange anomalie_massime.grib)
t2=$(grib_get -w count=2 -p stepRange anomalie_massime.grib)
NN1=$((t2-t1))

t3=$(grib_get -w count=1 -p stepRange anomalie_minime.grib)
t4=$(grib_get -w count=2 -p stepRange anomalie_minime.grib)
NN2=$((t4-t3))

printf -v NN1 "%02d" $NN1
printf -v NN2 "%02d" $NN2



# correction for vertical coordinates, dx and dy
# for the temperature anomalies the right indicatorOfParameter (25) has been specified; this allow GDAL to correctly read the temperature without apply a conversion from K to °C

grib_set -s indicatorOfParameter=25,deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 anomalie_massime.grib anomalieTemperatureMassime${gg}${mm}${yyyy}0000${NN1}${BB1}.grib

grib_set -s indicatorOfParameter=25,deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 anomalie_minime.grib anomalieTemperatureMinime${gg}${mm}${yyyy}0000${NN2}${BB2}.grib

grib_set -s deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 precacc_hires.grib precipitazioni${gg}${mm}${yyyy}0000${NN3}${BB3}.grib

grib_set -s deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 previ_tcorr_saturo_oggi.grib temperatura${gg}${mm}${yyyy}0000${NN4}${BB4}.grib

grib_set -s deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 wind_speed.grib vento${gg}${mm}${yyyy}0000${NN5}${BB5}.grib


# land-sea mask application

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       anomalieTemperatureMassime${gg}${mm}${yyyy}0000${NN1}${BB1}.grib anomalieTemperatureMassime_${gg}${mm}${yyyy}0000${NN1}${BB1}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       anomalieTemperatureMinime${gg}${mm}${yyyy}0000${NN2}${BB2}.grib anomalieTemperatureMinime_${gg}${mm}${yyyy}0000${NN2}${BB2}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       precipitazioni${gg}${mm}${yyyy}0000${NN3}${BB3}.grib precipitazioni_${gg}${mm}${yyyy}0000${NN3}${BB3}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       temperatura${gg}${mm}${yyyy}0000${NN4}${BB4}.grib temperatura_${gg}${mm}${yyyy}0000${NN4}${BB4}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       vento${gg}${mm}${yyyy}0000${NN5}${BB5}.grib vento_${gg}${mm}${yyyy}0000${NN5}${BB5}.grib


fi

#################################################################################
#run12

if [ "$hh" -ge 19 ]; then

b=12
c=36

for ((i=0, j=1; i<=11, j<=12; i+=1, j+=1)); do
    cdo max ${i}_previ.grib ${j}_previ.grib massime_oggi.grib
    if [ "$j" -ne "$b" ]; then
	grib_copy massime_oggi.grib ${j}_previ.grib
    fi
done

for ((i=24, j=25; i<=35, j<=36; i+=1, j+=1)); do
    cdo max ${i}_previ.grib ${j}_previ.grib massime_domani.grib
    if [ "$j" -ne "$c" ]; then
	grib_copy massime_domani.grib ${j}_previ.grib
    fi
done

rm -f *_previ.grib
grib_filter $SCRIPTS/filtra_step.txt previ_tcorr_saturo_oggi.grib

b=24
c=48

for ((i=12, j=13; i<=23, j<=24; i+=1, j+=1)); do
  
cdo min ${i}_previ.grib ${j}_previ.grib minime_oggi.grib

if [ "$j" -ne "$b" ]; then
grib_copy minime_oggi.grib ${j}_previ.grib
fi

done

for ((i=36, j=37; i<=47, j<=48; i+=1, j+=1)); do
  
cdo min ${i}_previ.grib ${j}_previ.grib minime_domani.grib

if [ "$j" -ne "$c" ]; then
grib_copy minime_domani.grib ${j}_previ.grib
fi
done



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

# number of bands BB?
BB1=$(grib_count anomalie_massime.grib)
BB2=$(grib_count anomalie_minime.grib)
printf -v BB1 "%02d" $BB1
printf -v BB2 "%02d" $BB2

t1=$(grib_get -w count=1 -p stepRange anomalie_massime.grib)
t2=$(grib_get -w count=2 -p stepRange anomalie_massime.grib)
NN1=$((t2-t1))


t3=$(grib_get -w count=1 -p stepRange anomalie_minime.grib)
t4=$(grib_get -w count=2 -p stepRange anomalie_minime.grib)
NN2=$((t4-t3))


printf -v NN1 "%02d" $NN1
printf -v NN2 "%02d" $NN2





# correction for vertical coordinates, dx and dy
# for the temperature anomalies the right indicatorOfParameter (25) has been specified; this allow GDAL to correctly read the temperature without apply a conversion from K to °C

grib_set -s indicatorOfParameter=25,deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 anomalie_massime.grib anomalieTemperatureMassime${gg}${mm}${yyyy}1200${NN1}${BB1}.grib
#anomalie_massime_finale.grib

grib_set -s indicatorOfParameter=25,deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 anomalie_minime.grib anomalieTemperatureMinime${gg}${mm}${yyyy}1200${NN2}${BB2}.grib
#anomalie_minime_finale.grib

grib_set -s deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 precacc_hires.grib precipitazioni${gg}${mm}${yyyy}1200${NN3}${BB3}.grib
#precacc_hires_finale.grib

grib_set -s deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 previ_tcorr_saturo_oggi.grib temperatura${gg}${mm}${yyyy}1200${NN4}${BB4}.grib
#previ_t_finale.grib

grib_set -s deletePV=1,earthIsOblate=1,iDirectionIncrement=17,jDirectionIncrement=17,resolutionAndComponentFlags=128 wind_speed.grib vento${gg}${mm}${yyyy}1200${NN5}${BB5}.grib
#wind_speed_finale.grib



# land-sea mask application

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       anomalieTemperatureMassime${gg}${mm}${yyyy}1200${NN1}${BB1}.grib anomalieTemperatureMassime_${gg}${mm}${yyyy}1200${NN1}${BB1}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       anomalieTemperatureMinime${gg}${mm}${yyyy}1200${NN2}${BB2}.grib anomalieTemperatureMinime_${gg}${mm}${yyyy}1200${NN2}${BB2}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       precipitazioni${gg}${mm}${yyyy}1200${NN3}${BB3}.grib precipitazioni_${gg}${mm}${yyyy}1200${NN3}${BB3}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       temperatura${gg}${mm}${yyyy}1200${NN4}${BB4}.grib temperatura_${gg}${mm}${yyyy}1200${NN4}${BB4}.grib

vg6d_transform --trans-type=metamorphosis --sub-type=maskvalid \
	       --maskbounds=0.1,1. --coord-file=mask_hires.grib --coord-format=grib_api \
	       vento${gg}${mm}${yyyy}1200${NN5}${BB5}.grib vento_${gg}${mm}${yyyy}1200${NN5}${BB5}.grib

fi

rm *_previ.grib
rm -f anomalieTemperatureMassime${gg}${mm}${yyyy}*.grib anomalieTemperatureMinime${gg}${mm}${yyyy}*.grib precipitazioni${gg}${mm}${yyyy}*.grib temperatura${gg}${mm}${yyyy}*.grib vento${gg}${mm}${yyyy}*.grib
