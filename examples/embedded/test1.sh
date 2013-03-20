#! /bin/bash

##
# test1.sh - Created by Timothy Morey on 3/12/2013
#
# This script contains the first test case for PISM's embedded mode, which is
# vaguely based on the examples/jako runs.
#


cp g5km_0_ftt.nc coarse-input.nc
cp jakofine_short.nc jakofine-input.nc


# Step 1: run pismr for 1 time step, writing out the state at both the beginning
# and end of the time step
mpiexec -n 8 pismr \
	-y 100 -step_count 1 \
	-i coarse-input.nc \
	-o coarse-output.nc \
	-step_record_file coarse-bc.nc -step_record_vars thk


# Step 2: run pismo until it catches up to pismr

# Step 2.1: Figure out the stop time of the coarse model in model years, and 
# store it in the variable STOP_TIME
STOP_TIME=`python extract_time.py coarse-output.nc`

mpiexec -n 8 pismo \
	-i jakofine-input.nc \
	-no_model_strip 10 \
	-ocean_kill jako.nc -cfbc -kill_icebergs -diffuse_bwat -thk_eff -sia_e 1.0 -ssa_sliding -topg_to_phi 5.0,30.0,-300.0,700.0 -plastic_pwfrac 0.98 -pseudo_plastic -pseudo_plastic_q 0.25 \
	-surface given,forcing -surface_given_file g5km_climate.nc -force_to_thk jako.nc \
	-ys 0 -ye $STOP_TIME \
	-skip -skip_max 10 \
	-coarse_grid_file coarse-bc.nc \
	-o jakofine-output.nc \
	-step_record_file jakofine-steps.nc -step_record_vars thk


# Step 3: combine pismr and pismo outputs to generate a new input for pismr

# TODO: repeat
