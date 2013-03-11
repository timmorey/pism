/*
 * CoarseGrid.hh
 *
 *  Created on: Mar 10, 2013
 *      Author: Timothy Morey
 */

#ifndef COARSEGRID_HH_
#define COARSEGRID_HH_

#include <petsc.h>
#include <string>
#include <vector>

class PIO;

class CoarseGrid
{
public:
  CoarseGrid(const std::string& filename);
  ~CoarseGrid();

public:
  PetscErrorCode SetAreaOfInterest(PetscReal xmin, PetscReal xmax,
                                   PetscReal ymin, PetscReal ymax);

protected:
  std::vector<PetscReal> x, y, z, t;
  int aoi_minxi, aoi_maxxi, aoi_minyi, aoi_maxyi;
  PIO* pio;

};

#endif /* COARSEGRID_HH_ */
