/*
 * CoarseGrid.hh
 *
 *  Created on: Mar 10, 2013
 *      Author: Timothy Morey
 */

#ifndef COARSEGRID_HH_
#define COARSEGRID_HH_

#include <map>
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

  PetscErrorCode Interpolate(const std::string& varname,
                             double x, double y, double z, double t,
                             double* value);

protected:
  std::vector<PetscReal> _X, _Y, _Z, _T;
  std::map<std::string, double*> _VarCache;
  int _AOIMinXi, _AOIMaxXi, _AOIMinYi, _AOIMaxYi;
  PIO* _Pio;

};

#endif /* COARSEGRID_HH_ */
