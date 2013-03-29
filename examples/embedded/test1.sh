#! /bin/bash

##
# test1.sh - Created by Timothy Morey on 3/12/2013
#
# This script contains the first test case for PISM's embedded mode, which is
# vaguely based on the examples/jako runs.
#


##
# Our input to the coarse grid is the fully spun-up 5km greenland grid, based on
# the examples/searise-greenland scripts.
cp g5km_0_ftt.nc coarse-input.nc

##
# Input for the fine grid is the fully spun-up 1km Jakobshavn grid based on the
# examples/jako scripts.  We ran the preprocess.sh, spinup.sh, and the first 
# bootstrapping run from century.sh.  This is a 1km grid with bed topography 
# based on the searise 1km dataset.
cp jakofine_short.nc jakofine-input.nc

##
# We're also currently using other data files produced by examples/jako scripts:
#   jako.nc - produced by quickjakosetup.sh
#   g5km_climate.nc - produced by preprocess.sh


# Step 1: run pismr for 1 time step, writing out the state at both the beginning
# and end of the time step
# The parameters for this run are based on the first control run given in 
# examples/searise-greenland/experiments.sh.
#mpiexec -n 8 pismr \
#  -bed_def lc -ssa_sliding -topg_to_phi 5.0,20.0,-300.0,700.0 -ocean_kill -climatic_mass_balance_cumulative \
#  -ocean constant -atmosphere searise_greenland -surface pdd,turn_into_anomaly -pdd_annualize \
#  -y 100 -step_count 1 \
#  -i coarse-input.nc \
#  -o coarse-output.nc \
#  -step_record_file coarse-bc.nc -step_record_vars thk,usurf,bmelt,vel_ssa,enthalpy


# Step 2: run pismo until it catches up to pismr

# Step 2.1: Figure out the stop time of the coarse model in model years, and 
# store it in the variable STOP_TIME
STOP_TIME=`python extract_time.py coarse-output.nc`

# Step 2.2: Run pismo until it catches up to pismr
# The parameters for this run are copied out of examples/jako/century.sh, with
# only minor modifications to fit it into the embedded framework.
mpiexec -n 8 pismo \
  -i jakofine-input.nc \
  -no_model_strip 10 \
  -ocean_kill jako.nc -cfbc -kill_icebergs -diffuse_bwat -thk_eff -sia_e 1.0 -ssa_sliding -topg_to_phi 5.0,30.0,-300.0,700.0 -plastic_pwfrac 0.98 -pseudo_plastic -pseudo_plastic_q 0.25 \
  -surface given,forcing -surface_given_file g5km_climate.nc -interp_ftt_thk coarse-bc.nc -ftt_mask_file jakofine-input.nc -force_to_thk_alpha 1.0 \
  -ys 0 -ye $STOP_TIME \
  -coarse_grid_file coarse-bc.nc \
  -o jakofine-output.nc \
  -step_record_file jakofine-steps.nc -step_record_vars thk,usurf,bmelt,vel_ssa

# Step 3: combine pismr and pismo outputs to generate a new input for pismr

# TODO: repeat
