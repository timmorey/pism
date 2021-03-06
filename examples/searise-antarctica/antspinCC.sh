#!/bin/bash

# Copyright (C) 2009-2011 Maria Martin and Ed Bueler
##################################################################################
# Spinup of Antarctic ice sheet model using data from Anne Le Brocq (from SeaRISE wiki).
# Uses PIK physics and enthalpy model (see Publications at www.pism-docs.org) 
# and modified configuration parameters with constant climate.  Uses constant
# and precip and a parameterization for artm as in Martin et al (2012).
# WARNING: especially for the last part, with finer resolutions, output is large!
##################################################################################


SCRIPTNAME="#(antspinCC.sh)"


echo "$SCRIPTNAME   run preprocess.sh before this..."
#inputs generated by preprocess.sh:
#  present day           = pism_Antarctica_5km.nc 					
#  glacial cycle forcing = pism_dT.nc pism_dSL.nc
#  future climate data   = ar4_ant_acab_anomaly_scalefactor_{1.0,1.5,2.0}.nc
#                          ar4_ant_artm_anomaly_scalefactor_{1.0,1.5,2.0}.nc

set -e  # exit on error

echo "$SCRIPTNAME   Constant-climate spinup-script using SeaRISE-Antarctica data"
echo "$SCRIPTNAME      and -ssa_sliding and -pik"
echo "$SCRIPTNAME   Run as './antspinCC.sh NN' for NN procs and 15km grid"


# naming files, directories, executables
RESDIR=
BOOTDIR=
PISM_EXEC=pismr
PISM_MPIDO="mpiexec -n "

# input data:
PISM_INDATANAME=${BOOTDIR}pism_Antarctica_5km.nc
PISM_TEMPSERIES=${BOOTDIR}pism_dT.nc
PISM_SLSERIES=${BOOTDIR}pism_dSL.nc

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "antspinup.sh 8" then NN = 8
  NN="$1"
fi
echo "$SCRIPTNAME              NN = $NN"
set -e  # exit on error

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi
DO=$PISM_DO

# grids
THIRTYKMGRID="-Mx 200 -My 200 -Lz 5000 -Lbz 2000 -Mz 41 -Mbz 16"
TWENTYKMGRID="-Mx 300 -My 300 -Lz 5000 -Lbz 2000 -Mz 41 -Mbz 16"
FIFTEENKMGRID="-Mx 400 -My 400 -Lz 5000 -Lbz 2000 -Mz 81 -Mbz 21"
TWELVEKMGRID="-Mx 500 -My 500 -Lz 5000 -Lbz 2000 -Mz 101 -Mbz 31"
TENKMGRID="-Mx 600 -My 600 -Lz 5000 -Lbz 2000 -Mz 101 -Mbz 31"
SEVENKMGRID="-Mx 900 -My 900 -Lz 5000 -Lbz 2000 -Mz 151 -Mbz 31"
FIVEKMGRID="-Mx 1200 -My 1200 -Lz 5000 -Lbz 2000 -Mz 201 -Mbz 51"

# skips:  
SKIPTHIRTYKM=10
SKIPTWENTYKM=10
SKIPFIFTEENKM=10
SKIPTWELVEKM=50
SKIPTENKM=100
SKIPSEVENKM=100
SKIPFIVEKM=200

SKIP=$SKIPFIFTEENKM

SIA_ENHANCEMENT="-sia_e 5.6"

#PIK-stuff; notes:
# 1)   '-pik' = '-cfbc -part_grid -part_redist -kill_icebergs'
# 2)   -meltfactor_pik 5e-3 is default when using -ocean pik
PIKPHYS="-ssa_method fd -ssa_e 0.6 -pik -eigen_calving -eigen_calving_K 2.0e18 -thickness_calving -calving_at_thickness 50.0"
#PIKPHYS_COUPLING="-atmosphere pik -ocean pik -meltfactor_pik 1.5e-2"
PIKPHYS_COUPLING="-atmosphere given -atmosphere_given_file $PISM_INDATANAME -surface simple -ocean pik -meltfactor_pik 1.5e-2"

# sliding related options:
PARAMS="-pseudo_plastic -pseudo_plastic_q 0.25 -plastic_pwfrac 0.97"
TILLPHI="-topg_to_phi 5.0,20.0,-300.0,700.0"
#TILLPHI="-topg_to_phi 5.0,20.0,-1000.0,0.0,10.0" # as in Martin et al 2012
FULLPHYS="-ssa_sliding -thk_eff $PARAMS $TILLPHI"


echo "$SCRIPTNAME             PISM = $PISM_EXEC"
echo "$SCRIPTNAME         FULLPHYS = $FULLPHYS"
echo "$SCRIPTNAME          PIKPHYS = $PIKPHYS"
echo "$SCRIPTNAME PIKPHYS_COUPLING = $PIKPHYS_COUPLING"


# #######################################
# bootstrap and SHORT smoothing run to 100 years
# #######################################
stage=smoothing
INNAME=$PISM_INDATANAME
RESNAME=${RESDIR}$stage.nc
RUNTIME=100 
echo
echo "$SCRIPTNAME  bootstrapping plus short SIA run for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -boot_file ${INNAME} $FIFTEENKMGRID \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING -ocean_kill \
	-y $RUNTIME -o $RESNAME"
$DO $cmd
#exit # <-- uncomment to stop here


# #######################################
# run with -no_mass (no surface change) on 15km for 200ka
# #######################################
stage=nomass
INNAME=$RESNAME
RESNAME=${RESDIR}$stage.nc
TSNAME=${RESDIR}ts_$stage.nc
RUNTIME=200000 
EXTRANAME=${RESDIR}extra_$stage.nc
expackage="-extra_times 0:1000:$RUNTIME -extra_vars bmelt,bwat,csurf,temppabase,diffusivity,hardav"
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME $PIKPHYS_COUPLING  \
	$SIA_ENHANCEMENT -no_mass \
  	-ys 0 -y $RUNTIME \
	-extra_file $EXTRANAME $expackage \
  	-o $RESNAME"
$DO $cmd
#exit # <-- uncomment to stop here


# #######################################
# run into steady state with constant climate forcing
# #######################################
stage=run
INNAME=$RESNAME
RESNAME=${RESDIR}$stage.nc
TSNAME=${RESDIR}ts_$stage.nc
RUNTIME=100000 
EXTRANAME=${RESDIR}extra_$stage.nc
exfilepackage="-extra_times 0:1000:$RUNTIME -extra_vars thk,usurf,cbase,cbar,mask,diffusivity,tauc,bmelt,bwat,temp "

echo
echo "$SCRIPTNAME  run into steady state with constant climate forcing for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -i $INNAME \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
	-ys 0 -y $RUNTIME \
	-ts_file $TSNAME -ts_times 0:1:$RUNTIME \
	-extra_file $EXTRANAME $exfilepackage \
	-o $RESNAME -o_size big"
$DO $cmd
#exit # <-- uncomment to stop here

# #######################################
## do a regridding to 10km:
# #######################################
#NN=64
#stage=run_regrid10km
#INNAME=${RESDIR}run.nc
#RESNAME=${RESDIR}$stage.nc
#OUTNAME=${OUTDIR}$stage.out
#TSNAME=${RESDIR}ts_$stage.nc
#RUNTIME=10000
#EXTRANAME=${RESDIR}extra_$stage.nc
#exfilepackage="-extra_times 0:500:$RUNTIME -extra_vars thk,usurf,cbase,cbar,mask,diffusivity,bmelt,bwat "

#echo
#echo "$SCRIPTNAME  run into steady state with constant climate forcing.. regridding to 10km"
#cmd="$PISM_MPIDO $NN $PISM_EXEC -skip $SKIPTENKM -regrid_file $INNAME -regrid_vars litho_temp,thk,enthalpy,bwat\
#	-boot_file $PISM_INDATANAME $TENKMGRID \
#	$SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
#	-ys 0 -y $RUNTIME \
#	-ts_file $TSNAME -ts_times 0:1:$RUNTIME \
#	-extra_file $EXTRANAME $exfilepackage \
#	-o $RESNAME -o_size big"
#	
#echo $DO $cmd 
#$DO $cmd >> $OUTNAME
##exit # <-- uncomment to stop here

# #######################################
## do a regridding to 6.7km:
# #######################################
#NN=128
#stage=run_regrid7km
#INNAME=${RESDIR}run.nc
#RESNAME=${RESDIR}$stage.nc
#OUTNAME=${OUTDIR}$stage.out
#TSNAME=${RESDIR}ts_$stage.nc
#RUNTIME=5000
#EXTRANAME=${RESDIR}extra_$stage.nc
#exfilepackage="-extra_times 0:500:$RUNTIME -extra_vars thk,usurf,cbase,cbar,mask,diffusivity,bmelt,bwat "

#echo
#echo "$SCRIPTNAME  run into steady state with constant climate forcing.. regridding to 6.7km"
#cmd="$PISM_MPIDO $NN $PISM_EXEC -skip $SKIPSEVENKM -regrid_file $INNAME -regrid_vars litho_temp,thk,enthalpy,bwat\
#	-boot_file $PISM_INDATANAME $SEVENKMGRID \
#	$SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
#	-ys 0 -y $RUNTIME \
#	-ts_file $TSNAME -ts_times 0:1:$RUNTIME \
#	-extra_file $EXTRANAME $exfilepackage \
#	-o $RESNAME -o_size big"
#	
#echo $DO $cmd 
#$DO $cmd >> $OUTNAME
##exit # <-- uncomment to stop here

# #######################################
## do a regridding to 5km:
# #######################################
#NN=128
#stage=run_regrid_5km
#INNAME=${RESDIR}run.nc
#RESNAME=${RESDIR}$stage.nc
#OUTNAME=${OUTDIR}$stage.out
#TSNAME=${RESDIR}ts_$stage.nc
#RUNTIME=1000
#EXTRANAME=${RESDIR}extra_$stage.nc
#exfilepackage="-extra_times 0:250:$RUNTIME -extra_vars thk,usurf,cbase,cbar,mask,diffusivity,bmelt,bwat "

#echo
#echo "$SCRIPTNAME  run into steady state with constant climate forcing"
#cmd="$PISM_MPIDO $NN $PISM_EXEC -skip $SKIPFIVEKM -regrid_file $INNAME -regrid_vars litho_temp,thk,enthalpy,bwat\
#	-boot_file $PISM_INDATANAME $FIVEKMGRID \
#	$SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
#	-ys 0 -y $RUNTIME \
#	-ts_file $TSNAME -ts_times 0:1:$RUNTIME \
#	-extra_file $EXTRANAME $exfilepackage \
#	-o $RESNAME -o_size big"
#	
#echo $DO $cmd 
#$DO $cmd >> $OUTNAME
##exit # <-- uncomment to stop here

echo
echo "$SCRIPTNAME  spinup done"


