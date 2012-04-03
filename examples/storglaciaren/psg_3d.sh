#!/bin/bash

# Copyright (C) 2011 Andy Aschwanden and Ed Bueler

set -e # exit on error

echo "# PISM Storglaciaren 3d Model"

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(psg_3d.sh)"
fi

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psg_3d.sh 8" then NN = 8
  NN="$1"
fi

echo "$SCRIPTNAME              NN = $NN"

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var DO is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=""    # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pismr -cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi

echo

PCONFIG=psg_config.nc

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -config_override $PCONFIG"


DATANAME=storglaciaren_3d.nc
PISM_DATANAME=pism_$DATANAME
INNAME=$PISM_DATANAME

COUPLER="-surface constant" # FIXME  should be using PSElevation as in flowline example

# 100 m grid
GRID="-Mx 38 -My 21 -Mz 51 -Mbz 1 -Lz 300 -z_spacing equal"
GS=100
SKIP=200

# 50 m grid
#GRID="-Mx 75 -My 41 -Mz 51 -Mbz 1 -Lz 300 -z_spacing equal"
#GS=50
#SKIP=200

# 20 m grid
#GRID="-Mx 186 -My 101 -Mz 301 -Mbz 1 -Lz 300 -z_spacing equal"
#GS=20
#SKIP=500

# 10 m grid
#GRID="-Mx 371 -My 201 -Mz 301 -Mbz 1 -Lz 300 -z_spacing equal"
#GS=10
#SKIP=500

REGRIDVARS="litho_temp,enthalpy,bwat,bmelt,thk"

RUNLENGTH=200
OUTNAME=psg_3d_${GS}m_${RUNLENGTH}a.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run for 1 a"
cmd="$PISM_MPIDO $NN $PISM -e 0.3 -skip -skip_max $SKIP -boot_file $INNAME $GRID \
  $COUPLER -ssa_sliding -plastic_phi 40 -thk_eff -y 100 -o $OUTNAME"
$PISM_DO $cmd


# 50 m grid
GRID="-Mx 75 -My 41 -Mz 51 -Mbz 1 -Lz 300 -z_spacing equal"
GS=50
SKIP=200


INNAME=$OUTNAME
RUNLENGTH=20
OUTNAME=psg_3d_${GS}m_${RUNLENGTH}a.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run for 1 a"
cmd="$PISM_MPIDO $NN $PISM -e 0.3 -skip -skip_max $SKIP -boot_file $INNAME $GRID \
  -regrid_file $INNAME -regrid_vars $REGRIDVARS $COUPLER -ssa_sliding -plastic_phi 40 -thk_eff -y 20 -o $OUTNAME"
$PISM_DO $cmd

# 20 m grid
GRID="-Mx 186 -My 101 -Mz 101 -Mbz 1 -Lz 300 -z_spacing equal"
GS=20
SKIP=500

INNAME=$OUTNAME
RUNLENGTH=5
OUTNAME=psg_3d_${GS}m_${RUNLENGTH}a.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run for 1 a"
cmd="$PISM_MPIDO $NN $PISM -e 0.3 -skip -skip_max $SKIP -boot_file $INNAME $GRID \
  -regrid_file $INNAME -regrid_vars $REGRIDVARS $COUPLER -ssa_sliding -plastic_phi 40 -thk_eff -y 20 -o $OUTNAME"
$PISM_DO $cmd

echo
echo "$SCRIPTNAME  done"

