#!/bin/bash

# Copyright (C) 2009-2012 The PISM Authors

# PISM SeaRISE Greenland spinup using modeled paleoclimate
#
# before using this script, run preprocess.sh to download and adjust metadata
# on SeaRISE "Present Day Greenland" master dataset
#
# recommended way to run with N processors is " ./spinup.sh N >& out.spinup & "
# which gives a viewable (with "less", for example) transcript in out.spinup
#

##
# This has been hacked by Timothy Morey, so that we can use this as an IO
# benchmark.  This required the following changes:
#   1) Reduce run lengths so that we focus on IO rather than computation, and
#      and also ensure that our runs are short enough to be practical on ranger.
#   2) For ranger friendliness, specify the number of nodes through the NN 
#      variable, rather than as a command line argument.  The reason for this is
#      that ranger doesn't require us to specify the number of processes for
#      each stage in this script as it is fixed for the batch.
#   3) Use consistent resolution for all tests, rather than starting at a coarse
#      resolution and progressing to a finer resolution for final runs.
#
# With these hacks, the usage of this script has changed.  No arguments are
# expected on the command line, and instead the behavior is controlled through
# the following environment variables:
#   PISM_MPIDO: The command used to launch the pism processes ("mpiexec -n <n>" 
#     (where <n> gives the number of processes) on most systems, but just "ibrun" 
#     on ranger).
#   PISM_RESOLUTION: The horizontal resolution for the simulation.  If set to 10, 
#     we will do a 10km resolution; if set to 5, we will do a 5km resolution; if
#     set to 1, we will do 1km resolution; for all other values, we will do a 20km 
#     resolution.
#   PISM_OFORMAT: The output format, either netcdf3, netcdf4_parallel, or pnetcdf.  
#     If not set, or given any other value, we will default to netcdf3.
#   PISM_DEFAULT_CONFIG: The path to the pism_config.nc file.  If this is not set, 
#     then we will not use the -config option, and PISM will try to find the file 
#     on its own.
##

# seconds per year, from UDUNITS
SECPERA=3.15569259747e7

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(spinup.sh)"
fi

echo
echo "# =================================================================================="
echo "# PISM SeaRISE Greenland: spinup"
echo "# =================================================================================="
echo

set -e  # exit on error

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n 8"
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

# preprocess.sh generates pism_*.nc files; run it first
if [ -n "${PISM_DATANAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_DATANAME = $PISM_DATANAME  (already set)"
else
  PISM_DATAVERSION=1.1
  PISM_DATANAME=pism_Greenland_5km_v$PISM_DATAVERSION.nc
fi
if [ -n "${PISM_TEMPSERIES:+1}" ] ; then
  echo "$SCRIPTNAME PISM_TEMPSERIES = $PISM_TEMPSERIES  (already set)"
else
  PISM_TEMPSERIES=pism_dT.nc
fi
if [ -n "${PISM_SLSERIES:+1}" ] ; then
  echo "$SCRIPTNAME   PISM_SLSERIES = $PISM_SLSERIES  (already set)"
else
  PISM_SLSERIES=pism_dSL.nc
fi
if [ -n "${PISM_CONFIG:+1}" ] ; then
  echo "$SCRIPTNAME     PISM_CONFIG = $PISM_CONFIG (already set)"
else
  PISM_CONFIG=searise_config.nc
fi

for INPUT in $PISM_DATANAME $PISM_TEMPSERIES $PISM_SLSERIES $PISM_CONFIG; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "$SCRIPTNAME           please run ./preprocess.sh, exiting"
    echo
    exit
  fi
done

INNAME=$PISM_DATANAME

# run lengths and starting time for paleo
SMOOTHRUNLENGTH=1 #100
NOMASSSIARUNLENGTH=1 #50000
PALEOSTARTYEAR=0 #-125000
FTTENDTIME=-1 #-100

# grids
VDIMS="-Lz 4000 -Lbz 2000"
COARSEVGRID="${VDIMS} -Mz 101 -Mbz 11 -z_spacing equal"
FINEVGRID="${VDIMS} -Mz 201 -Mbz 21 -z_spacing equal"
FINESTVGRID="${VDIMS} -Mz 401 -Mbz 41 -z_spacing equal"

TWENTYKMGRID="-Mx 76 -My 141 ${COARSEVGRID}"
TENKMGRID="-Mx 151 -My 281 ${FINEVGRID}"
FIVEKMGRID="-Mx 301 -My 561 ${FINESTVGRID}"
ONEKMGRID="-Mx 1501 -My 2801 ${FINIESTVGRID}"

# skips
SKIPTWENTYKM=10
SKIPTENKM=50
SKIPFIVEKM=200
SKIPONEKM=200

# defaults to coarse grid choices
COARSEGRID=$TWENTYKMGRID
FINEGRID=$TWENTYKMGRID
COARSESKIP=$SKIPTWENTYKM
FINESKIP=$SKIPTWENTYKM
CS=20 # km
FS=20 # km

COARSEENDTIME=0 #-5000 # BP

echo ""
if [ -z "$PISM_RESOLUTION" ] ; then
  PISM_RESOLUTION=20
fi
if [ $PISM_RESOLUTION -eq "10" ] ; then  # if user says "spinup.sh N 1" then MEDIUM:
  echo "$SCRIPTNAME grid: ALL RUNS ON 10km"
  echo "$SCRIPTNAME       WARNING: LARGE COMPUTATIONAL TIME"
  COARSEGRID=$TENKMGRID
  FINEGRID=$TENKMGRID
  COARSESKIP=$SKIPTENKM
  FINESKIP=$SKIPTENKM
  CS=10 # km
  FS=10 # km
elif [ $PISM_RESOLUTION -eq "5" ] ; then  # if user says "spinup.sh N 2" then FINE:
  echo "$SCRIPTNAME grid: ALL RUNS ON 5km"
  echo "$SCRIPTNAME       WARNING: VERY LARGE COMPUTATIONAL TIME"
  COARSEGRID=$FIVEKMGRID
  FINEGRID=$FIVEKMGRID
  COARSESKIP=$SKIPFIVEKM
  FINESKIP=$SKIPFIVEKM
  CS=5 # km
  FS=5 # km
elif [ $PISM_RESOLUTION -eq "1" ] ; then
  echo "$SCRIPTNAME grid: ALL RUNS ON 1km"
  echo "$SCRIPTNAME       WARNING: REALLY VERY SUPER LARGE COMPUTATION TIME"
  COARSEGRID=$ONEKMGRID
  FINEGRID=$ONEKMGRID
  COARSESKIP=$SKIPONEKM
  FINESKIP=$SKIPONEKM
  CS=1 # km
  FS=1 # km
else
    echo "$SCRIPTNAME grid: ALL RUNS ON 20km"
fi
echo ""

if [ "$PISM_OFORMAT" == "netcdf4_parallel" ] ; then
  echo "Using output format netcdf4_parallel"
elif [ "$PISM_OFORMAT" == "pnetcdf" ] ; then
  echo "Using output format pnetcdf"
elif [ "$PISM_OFORMAT" == "quilt" ] ; then
  echo "Using output format quilt"
else
  PISM_OFORMAT=netcdf3
  echo "Using output format netcdf3"
fi
echo ""

if [ -n "$PISM_DEFAULT_CONFIG" ] ; then
  PISM_DEFAULT_CONFIG="-config $PISM_DEFAULT_CONFIG"
fi

echo "$SCRIPTNAME     coarse grid = '$COARSEGRID' (= $CS km), with -skip = $COARSESKIP)"
echo "$SCRIPTNAME       fine grid = '$FINEGRID' (= $FS km), with -skip = $FINESKIP)"

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -config_override $PISM_CONFIG -climatic_mass_balance_cumulative -o_format $PISM_OFORMAT -log_summary $PISM_DEFAULT_CONFIG"

# coupler settings for pre-spinup
COUPLER_SIMPLE="-atmosphere searise_greenland -surface pdd -ocean_kill $INNAME"
# coupler settings for spin-up (i.e. with forcing)
COUPLER_FORCING="-atmosphere searise_greenland,delta_T -surface pdd -paleo_precip $PISM_TEMPSERIES -atmosphere_delta_T_file $PISM_TEMPSERIES -ocean constant,delta_SL -ocean_delta_SL_file $PISM_SLSERIES -ocean_kill $INNAME"
# coupler settings for spin-up (i.e. with forcing) and force-to-thickness
COUPLER_FTT="-atmosphere searise_greenland,delta_T -surface pdd,forcing -paleo_precip $PISM_TEMPSERIES -atmosphere_delta_T_file $PISM_TEMPSERIES -ocean constant,delta_SL -ocean_delta_SL_file $PISM_SLSERIES -ocean_kill $INNAME"

# default choices in parameter study; see Bueler & Brown (2009) re "tillphi"
TILLPHI="-topg_to_phi 5.0,20.0,-300.0,700.0"

FULLPHYS="-ssa_sliding -thk_eff ${TILLPHI}"

echo "$SCRIPTNAME      executable = '$PISM'"
echo "$SCRIPTNAME    full physics = '$FULLPHYS'"
echo "$SCRIPTNAME  simple coupler = '$COUPLER_SIMPLE'"
echo "$SCRIPTNAME forcing coupler = '$COUPLER_FORCING'"


# Run #1:
# bootstrap and do smoothing run to 100 years
PRE0NAME=g${CS}km_pre100.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run (for ${SMOOTHRUNLENGTH}a)"
cmd="$PISM_MPIDO $PISM -skip -skip_max  $COARSESKIP -boot_file $INNAME $COARSEGRID \
  $COUPLER_SIMPLE -y ${SMOOTHRUNLENGTH} -o $PRE0NAME"
$PISM_DO $cmd


# Run #2
# quick look at climate in 500 a recent period; see also delta_T in pism_dT.nc
#CLIMSTARTTIME=-500
#PRE0CLIMATE=g${CS}km_climate${CLIMSTARTTIME}a.nc
#PCLIM="${PISM_PREFIX}pclimate"
#echo
#echo "$SCRIPTNAME  running pclimate to show climate in modern period [${CLIMSTARTTIME} a,0 a], using current geometry and 10 year subintervals"
#cmd="$PISM_MPIDO $PCLIM -i $PRE0NAME $COUPLER_FORCING -times $CLIMSTARTTIME:10:0 -o $PRE0CLIMATE -o_format $PISM_OFORMAT -log_summary $PISM_DEFAULT_CONFIG"
#$PISM_DO $cmd


# Run #3
# run with -no_mass (no surface change) for 50ka
PRE1NAME=g${CS}km_steady.nc
EX1NAME=ex_${PRE1NAME}
EXTIMES=0:500:${NOMASSSIARUNLENGTH}
EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,bwp,csurf,hardav,mask" # check_stationarity.py can be applied to ex_${PRE1NAME}
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA run to achieve approximate temperature equilibrium, for ${NOMASSSIARUNLENGTH}a"
cmd="$PISM_MPIDO $PISM -skip -skip_max  $COARSESKIP -i $PRE0NAME $COUPLER_SIMPLE \
  -no_mass -ys 0 -y ${NOMASSSIARUNLENGTH} \
  -extra_file $EX1NAME -extra_vars $EXVARS -extra_times $EXTIMES -o $PRE1NAME"
$PISM_DO $cmd


# Run #4
# pre-spinup done; ready to use paleoclimate forcing for real spinup ...

EXSTEP=500
EXFSTEP=10
TSSTEP=yearly

STARTTIME=0 #$PALEOSTARTYEAR

ENDTIME=1 #$COARSEENDTIME

ET=$(($ENDTIME/-1))
OUTNAME=g${CS}km_m${ET}a.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,bwp,csurf,hardav,mask,dHdt,cbase,tauc,thk,topg,usurf,climatic_mass_balance_cumulative"
echo
echo "$SCRIPTNAME  paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation, from $PALEOSTARTYEAR a to ${ENDTIME}a"
cmd="$PISM_MPIDO $PISM -skip -skip_max  $COARSESKIP -i $PRE1NAME $FULLPHYS -bed_def lc $COUPLER_FORCING \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

# exit    # uncomment to stop here


# Run #5
# ######################################
# "regular" run
# ######################################

STARTTIME=0 #$ENDTIME
ENDTIME=1 # BP
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_0.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
echo
echo "$SCRIPTNAME  regular run"
echo "$SCRIPTNAME  regrid to fine grid and do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation,"
echo "$SCRIPTNAME      from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $PISM -skip -skip_max  $FINESKIP -boot_file $INNAME $FINEGRID $FULLPHYS \
     -bed_def lc $COUPLER_FORCING \
     -regrid_file $STARTNAME -regrid_vars litho_temp,thk,enthalpy,bwat,bmelt -regrid_bed_special  \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd


# Run #6
# ######################################
# "force-to-thickness" run
# ######################################

STARTTIME=-2 #$COARSEENDTIME
ENDTIME=$FTTENDTIME
ET=$(($FTTENDTIME/-1))
OUTNAME=g${FS}km_m${ET}a_ftt.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$FTTENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$FTTENDTIME
echo
echo "$SCRIPTNAME  force-to-thickness run"
echo "$SCRIPTNAME  regrid to fine grid and do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation and modified surface mass balance,"
echo "$SCRIPTNAME      from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $PISM -skip -skip_max  $FINESKIP -boot_file $INNAME $FINEGRID $FULLPHYS \
     -bed_def lc $COUPLER_FTT \
     -force_to_thk $INNAME -force_to_thk_alpha 0.005 \
     -regrid_file $STARTNAME -regrid_vars litho_temp,thk,enthalpy,bwat -regrid_bed_special  \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $FTTENDTIME -o $OUTNAME"
$PISM_DO $cmd


# Run #7

STARTTIME=$FTTENDTIME
ENDTIME=0
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_0_ftt.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXFSTEP:$ENDTIME
echo
echo "$SCRIPTNAME  force-to-thickness run finishes with frequent diagnostic information"
echo "$SCRIPTNAME  do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation and modified surface mass balance,"
echo "$SCRIPTNAME      from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $PISM -skip -skip_max  $FINESKIP -i $STARTNAME $FULLPHYS \
     -bed_def lc $COUPLER_FTT \
     -force_to_thk $INNAME -force_to_thk_alpha 0.005 \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd


# With hacks to reduce run lengths, we no longer write out any extra data.  As 
# such, these postprocessing runs have nothing to do and will only produce 
# errors.

#echo
#echo "$SCRIPTNAME  some postprocessing"
#echo
## calculate yearly-averages of acab and dHdt using ncap2 sleight of hand.
#ncap2_script="*sz_idt=time.size(); acab[\$time,\$x,\$y]= 0.f; dHdt[\$time,\$x,\$y]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {acab(idt,:,:)=(climatic_mass_balance_cumulative(idt,:,:)-climatic_mass_balance_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA; dHdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA;}"
#
#$PISM_DO ncap2 -O -s "$ncap2_script" $EXNAME $EXNAME
#
#echo
## adjust meta data for new fields
#$PISM_DO ncatted -a units,acab,o,c,"m year-1" -a units,dHdt,o,c,"m year-1" \
#      -a long_name,acab,o,c,"surface mass balance" \
#      -a long_name,dHdt,o,c,"rate of change of ice thickness" \
#      -a grid_mapping,acab,o,c,"mapping" \
#      -a grid_mapping,dHdt,o,c,"mapping" \
#      -a cell_methods,acab,o,c,"time: mean (interval: $EXFSTEP years)" \
#      -a cell_methods,dHdt,o,c,"time: mean (interval: $EXFSTEP years)" $EXNAME

#echo
## now extract last acab record
#cmd="ncks -A -v acab -d time,$ENDTIME. $EXNAME $OUTNAME"
#$PISM_DO $cmd

echo
echo "$SCRIPTNAME  spinup done"

