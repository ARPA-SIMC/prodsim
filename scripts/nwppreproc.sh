#!/bin/bash

# script for preliminary preprocessing of NWP data

# process on hourly intervals
DT='0 01'

# recompute statistically processed values on intervals of equal
# length $DT
# accumulated fields
vg6d_transform --comp-stat-proc=0:0 --comp-step=$DT $1 $1.acc
# averaged fields
vg6d_transform --comp-stat-proc=0:0 --comp-step=$DT $1 $1.avg

# compute derived quantities (see variable table in
# /usr/share/wreport/dballe.txt)
# wind speed and direction, relative humidity
vg6d_transform --output-variable-list=B11001,B11002,B13003 $1 $1.newvars

