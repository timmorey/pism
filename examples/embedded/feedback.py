##
# feedback.py - Created by Timothy Morey on 3/30/2013
#
# This script extracts relevant information from the output of the embedded
# model and interpolates for feedback into the coarse model.  We currently
# employ some simplifying assumptions:
#   - Both files contain a single time slice of data and they both contain the
#     same time slice, so we don't perform any temporal interpolation.
#   - Both files use the same CRS.
#


import netCDF4
import numpy
import os
import sys


def feedback(finefile, coarsefile):
    ncfine = netCDF4.Dataset(finefile, 'r')
    nccoarse = netCDF4.Dataset(coarsefile, 'r+')

    fx = ncfine.variables['x'][:]
    fy = ncfine.variables['y'][:]
    fz = ncfine.variables['z'][:]
    fthk = ncfine.variables['thk'][:]

    cx = nccoarse.variables['x'][:]
    cy = nccoarse.variables['y'][:]
    cz = nccoarse.variables['z'][:]
    cthk = nccoarse.variables['thk'][:]

    for i in range(len(cx)):
        x = cx[i]
        for j in range(len(cy)):
            y = cy[j]
            if x >= fx[0] and x <= fx[-1] and y >= fy[0] and y <= fy[-1]:
                xi1, xi2, x1, x2 = None, None, None, None
                yi1, yi2, y1, y2 = None, None, None, None

                for xi in range(1, len(fx)):
                    if fx[xi-1] <= x and x <= fx[xi]:
                        xi1, xi2 = xi-1, xi
                        x1, x2 = fx[xi1], fx[xi2]
                        break

                for yi in range(1, len(fy)):
                    if fy[yi-1] <= y and y <= fy[yi]:
                        yi1, yi2 = yi-1, yi
                        y1, y2 = fy[yi1], fy[yi2]
                        break

                values = [fthk[0, xi1, yi1],
                          fthk[0, xi1, yi2],
                          fthk[0, xi2, yi1],
                          fthk[0, xi2, yi2]]
                weights = [((x2-x)*(y2-y)/(x2-x1)*(y2-y1)),
                           ((x2-x)*(y-y1)/(x2-x1)*(y2-y1)),
                           ((x-x1)*(y2-y)/(x2-x1)*(y2-y1)),
                           ((x-x1)*(y-y1)/(x2-x1)*(y2-y1))]

                value = 0.0
                for k in range(4):
                    value += values[k] * weights[k]

                cthk[0, i, j] = value

    ncfine.close()

    nccoarse.variables['thk'][:] = cthk
    nccoarse.close()
                

if __name__ == '__main__':
    coarsefile = None
    finefile = None

    if len(sys.argv) > 2:
        finefile = sys.argv[1]
        coarsefile = sys.argv[2]

    if (coarsefile is not None and os.path.exists(coarsefile) and
        finefile is not None and os.path.exists(finefile)):
        feedback(finefile, coarsefile)
    else:
        print 'Error: Invalid input'
