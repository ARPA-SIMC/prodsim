#!/bin/bash

ver_uvrtg() {

    grib_get -p uvRelativeToGrid $1 | sort -u

    
}

interp_er() {

    vg6d_transform --trans-type=inter --sub-type=near --type=regular_ll \
		   --x-min=9. --y-min=43.5 --x-max=14. --y-max=45.5 \
		   --nx=251 --ny=101 \
		   $1 $2
    # --output-variable-list=B11002
    }


prepare_const_g1() {

    # vertical coordinate from half to full levels
    vg6d_transform --trans-type=vertint --sub-type=linear \
		   --trans-level-type=105,,105,105 --component-flag=1 \
		   in/hhl_i2.grib hfl_i2.grib

    # ver_uvrtg hfl_i2.grib

    # extract orography from last half level, change 66 to something else
    # if needed!
    grib_copy -w level=66 in/hhl_i2.grib zs.grib
    grib_set -s typeOfLevel=surface zs.grib zsurf.grib
    # add it to vertical full levels, needed for interpolation to height
    # from ground
    cat zsurf.grib >> hfl_i2.grib
    # interpolate to E-R area
    # interp_er hfl_i2.grib hfl_i2_er.grib
    rm -f zs.grib zsurf.grib

    }


prepare_const_g2() {

    grib_set -s typeOfSecondFixedSurface=255 in/hhl_i2.grib hhl_i2.tmp.grib
    
    # vertical coordinate from half to full levels
    vg6d_transform --trans-type=vertint --sub-type=linear \
		   --trans-level-type=150,,150,150 --component-flag=1 \
		   hhl_i2.tmp.grib hfl_i2.tmp.grib

    # this is required if data is grib1 (level=105)
    grib_set -s typeOfFirstFixedSurface=105,typeOfSecondFixedSurface=105 \
	     hfl_i2.tmp.grib hfl_i2.grib
    # ver_uvrtg hfl_i2.grib

    # extract orography from last half level, change 66 to something else
    # if needed!
    grib_copy -w level=66 hhl_i2.tmp.grib zs.grib
    rm -f hhl_i2.tmp.grib hfl_i2.tmp.grib hfl_i2.tmp.grib
    grib_set -s typeOfFirstFixedSurface=1,typeOfSecondFixedSurface=255 \
	     zs.grib zsurf.grib
    # add it to vertical full levels, needed for interpolation to height
    # from ground
    cat zsurf.grib >> hfl_i2.grib
    # interpolate to E-R area
    # interp_er hfl_i2.grib hfl_i2_er.grib
    rm -f zs.grib zsurf.grib

    }

LOG4C_PRIORITY=warn

set -x
rm *.grib
# https://arpa-simc.github.io/dballe/general_ref/ltypes.html

#### to be done once in the life ####

prepare_const_g2

#### to be done for every time level ####

# destagger data
vg6d_transform --a-grid in/uv_i2.grib uvdestag_i2.grib

# ver_uvrtg uvdestag_i2.grib


# interpolate to levels above mean sea (not what desired!)
#vg6d_transform --trans-type=vertint --sub-type=linear \
#	       --trans-level-type=105,105,102, \
#	       --trans-level-list=100000,150000,200000,250000,300000 \
#	       --coord-file=hfl_i2.grib --coord-format=grib_api \
#	       uvdestag_i2.grib uvzabs_i2.grib


# ver_uvrtg uvzabs_i2.grib

# vertical interpolation to levels (in mm) over surface,
# --anavariable-list=B10007 not needed since height of orography is
# provided in coord-file and not in data
vg6d_transform --trans-type=vertint --sub-type=linear \
	       --trans-level-type=105,105,103,  \
	       --trans-level-list=100000,150000,200000,250000,300000 \
	       --coord-file=hfl_i2.grib --coord-format=grib_api \
	       uvdestag_i2.grib uvzrel_i2.grib


# interpolate to E-R area and compute wind speed
interp_er "--output-variable-list=B11002 uvzrel_i2.grib" uvzrel_i2_er.grib

# ver_uvrtg uvzrel_i2.grib

