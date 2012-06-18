#!/bin/bash

./sg_create_flowline.py
./temp_forcing.py --e 2 --te 50 --n 50 delta_T.nc
DATANAME=storglaciaren_flowline.nc
PISM_DATANAME=pism_$DATANAME
flowline.py -e -o $PISM_DATANAME $DATANAME
# config file
CDLCONFIG=psg_config.cdl
PCONFIG=psg_config.nc
ncgen -o $PCONFIG $CDLCONFIG
