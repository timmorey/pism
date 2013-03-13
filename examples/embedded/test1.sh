#! /bin/bash

##
# test1.sh - Created by Timothy Morey on 3/12/2013
#
# This script contains the first test case for PISM's embedded mode, which is
# vaguely based on the examples/jako runs.
#

cp g5km_0_ftt.nc coarse-input.nc

# Step 1: run pismr for 1 time step, writing out the state at both the beginning
# and end of the time step
mpiexec -n 8 pismr \
	-i coarse-input.nc \
	-y 100 -step_count 1 \
	-o_size none \
	-save_file coarse-output.nc -save_times 0:0.01:100 -save_size medium

# Step 2: run pismo until it catches up to pismr

# Step 2.1: Figure out the stop time of the coarse model in model years, and 
# store it in the variable STOP_TIME
STOP_TIME=`python extract_time.py coarse-output.nc`

# TODO: run pismo until it catches up where pismr stopped
# TODO: use "-coarse_grid_file coarse-output.nc" to use coarse model state for boundary conditions

# TODO: combine pismr and pismo outputs to generate a new input for pismr

# TODO: repeat
