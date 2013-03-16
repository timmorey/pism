##
# extract_time.py - Created by Timothy Morey on 3/13/2013
#
# This script performs a single simple operation: it opens a netcdf file, looks
# for a 'time' variable, and prints the greatest time value to stdout.
#

import netCDF4
import numpy
import os
import sys

def findGreatestTime(ncfile):
    maxt = -sys.float_info.max
    nc = netCDF4.Dataset(ncfile, 'r')
    tvar = nc.variables['time']
    for t in tvar:
        maxt = max(maxt, t)
    return maxt

def secondsToYears(t):
    return t / 60 / 60 / 24 / 365

if __name__ == '__main__':
    infile = None

    if len(sys.argv) > 1:
        infile = sys.argv[1]

    if infile is not None and os.path.exists(infile):
        print secondsToYears(findGreatestTime(infile))
    else:
        print 'Error: Invalid input'
