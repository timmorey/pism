#!/usr/bin/env python

## Copyright (C) 2011 The PISM Authors

## script to generate figure: results from SeaRISE experiments
## usage:  if UAFX_G_D3_C?_??.nc are result NetCDF files then do
##   $ slr_show.py -m UAFX

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

from numpy import zeros
import pylab as plt
from optparse import OptionParser

parser = OptionParser()
parser.usage = "usage: %prog [options]"
parser.description = "A script for PISM output files to show time series plots using pylab."
parser.add_option("-m", "--model",dest="model",
                  help="choose experiment, default UAF1",default="UAF1")


(options, args) = parser.parse_args()
model = options.model

# first name in this list is CONTROL
NCNAMES = [model + "_G_D3_C1_E0.nc",model + "_G_D3_C2_E0.nc",model + "_G_D3_C3_E0.nc",model + "_G_D3_C4_E0.nc",model + "_G_D3_C1_S1.nc",model + "_G_D3_C1_S2.nc",model + "_G_D3_C1_S3.nc",model + "_G_D3_C1_M1.nc",model + "_G_D3_C1_M2.nc",model + "_G_D3_C1_M3.nc"]

# labels
labels = ["AR4 A1B","AR4 A1B 1.5x","AR4 A1B 2x","2x basal sliding","2.5x basal sliding","3x basal sliding","2 m/a bmr","20 m/a bmr","200 m/a bmr"]

print "control run name is " + NCNAMES[0]

n = len(NCNAMES)
nc0 = CDF(NCNAMES[0], 'r')
try:
  t = nc0.variables['tseries'][:]
except:
  t = nc0.variables['time'][:]
nc0.close()

ivol = zeros((len(t),n))
ivolshift = zeros((len(t),n-1))

for j in range(n):
  nc = CDF(NCNAMES[j], 'r')
  ivol[:,j] = nc.variables['ivol'][:]
  nc.close()

for j in range(n-1):
  ivolshift[:,j] = ivol[:,j+1] - ivol[:,0]

# "2,850,000 km3 of ice were to melt, global sea levels would rise 7.2 m"
scale = 7.2 / 2.850e6
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t,-(ivolshift/1.0e9)*scale)
ax.set_xlabel('years from 2004')
ax.set_ylabel('sea level rise relative to control (m)')
ax.legend(labels,loc='upper left')
ax.grid(True)

plt.show()
