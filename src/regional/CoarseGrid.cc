/*
 * CoarseGrid.cc
 *
 *  Created on: Mar 10, 2013
 *      Author: Timothy Morey
 */

#include "CoarseGrid.hh"
#include "PIO.hh"

#include <petsc.h>
#include <string>
#include <vector>


CoarseGrid::CoarseGrid(const std::string& filename)
  : aoi_minxi(0),
    aoi_maxxi(0),
    aoi_minyi(0),
    aoi_maxyi(0),
    pio(0)
{
  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  unsigned int xlen, ylen, zlen, tlen;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  pio = new PIO(comm, rank, "guess_mode");
  if(0 != pio->open(filename, PISM_NOWRITE)) {
    fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to open '%s'.\n", filename.c_str());
    delete pio;
    pio = 0;
  }

  if(pio) {
    if(0 != pio->inq_dimlen("x", xlen) ||
       0 != pio->inq_dimlen("y", ylen) ||
       0 != pio->inq_dimlen("z", zlen) ||
       0 != pio->inq_dimlen("t", tlen)) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to retrieve dimlen.\n");
      pio->close();
      delete pio;
      pio = 0;
    }
  }

  if(pio) {
    x.reserve(xlen);
    y.reserve(ylen);
    z.reserve(zlen);
    t.reserve(tlen);

    aoi_maxxi = xlen - 1;
    aoi_maxyi = ylen - 1;

    if(0 != pio->get_1d_var("x", 0, xlen, x) ||
       0 != pio->get_1d_var("y", 0, ylen, y) ||
       0 != pio->get_1d_var("z", 0, zlen, z) ||
       0 != pio->get_1d_var("t", 0, tlen, t)) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to read dimension scales.\n");
      pio->close();
      delete pio;
      pio = 0;
    }
  }
}

PetscErrorCode CoarseGrid::SetAreaOfInterest(PetscReal xmin, PetscReal xmax,
                                             PetscReal ymin, PetscReal ymax)
{
  PetscErrorCode retval = 0;

  // TODO: We assume that the region of interest is fully contained in the
  // bounds of this coarse grid.

  for(int i = 1; i < x.size(); i++) {
    if(xmin < x[i]) {
      aoi_minxi = i - 1;
      break;
    }
  }

  for(int i = x.size() - 2; i >= 0; i--) {
    if(xmax > x[i]) {
      aoi_maxxi = i + 1;
      break;
    }
  }

  for(int i = 1; i < y.size(); i++) {
    if(ymin < y[i]) {
      aoi_minyi = i - 1;
      break;
    }
  }

  for(int i = y.size() - 2; i >= 0; i--) {
    if(ymax > y[i]) {
      aoi_maxyi = i + 1;
      break;
    }
  }

  return retval;
}
