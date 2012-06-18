#!/usr/bin/env python

# Copyright (C) 2011 Andy Aschwanden

try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
NC = netCDF.Dataset

import numpy as np
from optparse import OptionParser

# Set up the option parser
parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "Script adds ocean forcing to HIRHAM atmosphere/surface forcing file. Sets a constant, spatially-uniform basal melt rate of b_a before time t_a, and b_e after time t_a."
parser.add_option("--a",dest="Ta",
                  help="start temp",default=0)
parser.add_option("--e",dest="Te",
                  help="end temp",default=2)
parser.add_option("--ta",dest="ta",
                  help="start year",default=0)
parser.add_option("--te",dest="te",
                  help="end year",default=100)
parser.add_option("--n",dest="n",
                  help="no steps",default=100)

(options, args) = parser.parse_args()
Ta = float(options.Ta)
Te = float(options.Te)
ta = float(options.ta)
te = float(options.te)
n = float(options.n)

infile = args[0]

nc = NC(infile,'w')

nc.createDimension("time", size=n+1)
time_var = nc.createVariable("time", 'f', dimensions=("time",))
time_var.units = "years";
time_var[:] = np.linspace(ta,te,n+1)

delta_T = np.linspace(Ta,Te,n+1)

fill_value = -2e33

temp_var = nc.createVariable("delta_T", 'f', dimensions=("time",))
temp_var.units = "K";
temp_var.longname = "Temperature (variation from present)"
temp_var.standard_name = "land_ice_temperature_at_firn_base"
temp_var[:] = delta_T

nc.close()
