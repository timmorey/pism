#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# see README for role of this script

# This script uses NCO (http://nco.sourceforge.net/).

set -e  # exit on error

echo
echo "# PDDTUNE example: preprocessing"
echo

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATANAME=Greenland_5km_v1.1.nc
PISMDATA=pism_$DATANAME
wget -nc http://websrv.cs.umt.edu/isis/images/a/a5/$DATANAME

ncks -O $DATANAME $PISMDATA  # copies over, but with history of this copy

# adjust metadata
ncrename -O -v x1,x -v y1,y -d x1,x -d y1,y $PISMDATA
ncrename -O -v time,t -d time,t $PISMDATA
ncrename -O -v usrf,usurf $PISMDATA

# convert from water-equiv to ice-equiv thickness change rates
#   using ice density 910.0 kg m-3
ncap2 -O -s "precip=presprcp*(1000.0/910.0)" $PISMDATA $PISMDATA
ncks -O -x -v presprcp $PISMDATA $PISMDATA  # delete old field, for clarity
ncatted -O -a units,precip,m,c,"m a-1" $PISMDATA
ncatted -O -a long_name,precip,m,c,"ice-equivalent precipitation rate" $PISMDATA
ncatted -a standard_name,precip,d,, $PISMDATA
ncap2 -O -s "newrunoff=runoff*(1000.0/910.0)" $PISMDATA $PISMDATA
ncks -O -x -v runoff $PISMDATA $PISMDATA  # delete old field, for clarity
ncrename -O -v newrunoff,runoff $PISMDATA
ncatted -O -a units,runoff,m,c,"m a-1" $PISMDATA
ncatted -O -a long_name,runoff,m,c,"ice-equivalent rate of mass loss through runoff; believed to include rain" $PISMDATA
ncatted -a standard_name,runoff,d,, $PISMDATA
ncap2 -O -s "newsmb=smb*(1000.0/910.0)" $PISMDATA $PISMDATA
ncks -O -x -v smb $PISMDATA $PISMDATA  # delete old field, for clarity
ncrename -O -v newsmb,smb $PISMDATA
ncatted -O -a units,smb,m,c,"m a-1" $PISMDATA
ncatted -O -a long_name,smb,m,c,"ice-equivalent surface mass balance rate" $PISMDATA
ncatted -a standard_name,smb,d,, $PISMDATA

# delete some more fields which we will not use, to avoid confusion
ncks -O -x -v airtemp2m,surftemp $PISMDATA $PISMDATA  # we will use Fausto parameterization

# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMDATA
echo "  PISM-readable (-boot_file) file $PISMDATA created from $DATANAME"

# form NetCDF version of basic config file
CONFIG=base_config
echo
echo "  generating parameter-override (-config_override) $CONFIG.nc file from $CONFIG.cdl"
rm -rf $CONFIG.nc   # remove old if present
ncgen -o $CONFIG.nc $CONFIG.cdl
