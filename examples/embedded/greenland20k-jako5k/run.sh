#! /bin/bash

##
# run.sh - Created by Timothy Morey on 4/1/2013
#
# This script aims to produce a fast-moving embedded test case for pism, to
# help with debugging.
#i

cp g20km_0_ftt.nc coarse-input.nc
cp jako.nc jako-input.nc

# Don't append to an old step file.  If the step files exist, move them to
# force new ones to be created.
if [ -f coarse-steps.nc ];
then
  mv coarse-steps.nc coarse-steps.nc~
fi

if [ -f jako-steps.nc ];
then
  mv jako-steps.nc jako-steps.nc~
fi


SIM_START=0
SIM_STOP=1

START_TIME=$SIM_START

while [ 1 -eq `echo "$START_TIME < $SIM_STOP" | bc` ]
do

  # Step 1: run pismr for 1 time step, writing out the state at both the beginning
  # and end of the time step
  # The parameters for this run are based on the first control run given in 
  # examples/searise-greenland/experiments.sh.
  mpiexec -n 8 pismr \
    -ssa_sliding \
    -ocean constant -ocean_kill -atmosphere searise_greenland -surface pdd,turn_into_anomaly -pdd_annualize \
    -ys $START_TIME -ye $SIM_STOP -step_count 1 \
    -i coarse-input.nc \
    -o coarse-output.nc \
    -step_record_file coarse-bc.nc -step_record_vars thk,usurf,bmelt,vel_ssa,enthalpy 

  # Step 2: run pismo until it catches up to pismr

  # Step 2.1: Figure out the stop time of the coarse model in model years, and 
  # store it in the variable STOP_TIME
  STOP_TIME=`python ../extract_time.py coarse-output.nc`

  # Step 2.2: Run pismo until it catches up to pismr
  # The parameters for this run are copied out of examples/jako/century.sh, with
  # only minor modifications to fit it into the embedded framework.
  mpiexec -n 8 pismo \
    -i jako-input.nc \
    -no_model_strip 4 \
    -ssa_sliding \
    -surface pdd,forcing -interp_ftt_thk coarse-bc.nc -ftt_mask_file jako-input.nc -force_to_thk_alpha 1.0 \
    -ys $START_TIME -ye $STOP_TIME \
    -coarse_grid_file coarse-bc.nc \
    -o jako-output.nc \
    -step_record_file jako-steps.nc -step_record_vars thk,usurf,bmelt,vel_ssa -append_steps

  # Step 3: combine pismr and pismo outputs to generate a new input for pismr
  python ../feedback.py jako-output.nc coarse-output.nc

  START_TIME=$STOP_TIME
  cp coarse-output.nc coarse-input.nc
  cp jako-output.nc jako-input.nc

done
