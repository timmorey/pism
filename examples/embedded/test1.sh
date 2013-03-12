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
mpiexc -n 8 pismr \
	-i coarse-input.nc \
	-y 10 -step_count 0 \
	-o_size none \
	-save_file coarse-output.nc -save_times 0:0.1:10 -save_size medium

# Step 2: run pismo until it catches up to pismr
# TODO: figure out the time that pismr ran to, so that we can stop pismo at the same place
# TODO: run pismo until it catches up where pismr stopped
# TODO: use "-coarse_grid_file coarse-output.nc" to use coarse model state for boundary conditions

# TODO: combine pismr and pismo outputs to generate a new input for pismr

# TODO: repeat
