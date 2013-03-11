/*
 * CoarseGrid.cc
 *
 *  Created on: Mar 10, 2013
 *      Author: Timothy Morey
 */

#include "CoarseGrid.hh"
#include "PIO.hh"

#include <map>
#include <petsc.h>
#include <stdlib.h>
#include <string>
#include <vector>


CoarseGrid::CoarseGrid(const std::string& filename)
  : _AOIMinXi(0),
    _AOIMaxXi(0),
    _AOIMinYi(0),
    _AOIMaxYi(0),
    _Pio(0) {

  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  unsigned int xlen, ylen, zlen, tlen;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  printf("Opening '%s'...\n", filename.c_str());
  _Pio = new PIO(comm, rank, "netcdf3");
  if(0 != _Pio->open(filename, PISM_NOWRITE)) {
    fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to open '%s'.\n", filename.c_str());
    delete _Pio;
    _Pio = 0;
  }

  printf("Reading dimension sizes from coarse grid file...\n");
  if(_Pio) {
    if(0 != _Pio->inq_dimlen("x", xlen) ||
       0 != _Pio->inq_dimlen("y", ylen) ||
       0 != _Pio->inq_dimlen("z", zlen) ||
       0 != _Pio->inq_dimlen("time", tlen)) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to retrieve dimlen.\n");
      _Pio->close();
      delete _Pio;
      _Pio = 0;
    }
  }

  printf("Reading coordinate values from coarse grid file...\n");
  if(_Pio) {
    _X.reserve(xlen);
    _Y.reserve(ylen);
    _Z.reserve(zlen);
    _T.reserve(tlen);

    _AOIMaxXi = xlen - 1;
    _AOIMaxYi = ylen - 1;

    if(0 != _Pio->get_1d_var("x", 0, xlen, _X) ||
       0 != _Pio->get_1d_var("y", 0, ylen, _Y) ||
       0 != _Pio->get_1d_var("z", 0, zlen, _Z) ||
       0 != _Pio->get_1d_var("time", 0, tlen, _T)) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to read dimension scales.\n");
      _Pio->close();
      delete _Pio;
      _Pio = 0;
    }
  }

  printf("Finished initializing CoarseGrid.\n");
}

CoarseGrid::~CoarseGrid() {
  if(_Pio) {
    _Pio->close();
    delete _Pio;
    _Pio = 0;
  }

  map<std::string, double*>::iterator mapiter;
  for(mapiter = _VarCache.begin(); mapiter != _VarCache.end(); mapiter++) {
    delete mapiter->second;
  }
}

PetscErrorCode CoarseGrid::SetAreaOfInterest(PetscReal xmin, PetscReal xmax,
                                             PetscReal ymin, PetscReal ymax) {
  PetscErrorCode retval = 0;

  // TODO: We assume that the region of interest is fully contained in the
  // bounds of this coarse grid.

  for(size_t i = 1; i < _X.size(); i++) {
    if(xmin < _X[i]) {
      _AOIMinXi = i - 1;
      break;
    }
  }

  for(int i = _X.size() - 2; i >= 0; i--) {
    if(xmax > _X[i]) {
      _AOIMaxXi = i + 1;
      break;
    }
  }

  for(size_t i = 1; i < _Y.size(); i++) {
    if(ymin < _Y[i]) {
      _AOIMinYi = i - 1;
      break;
    }
  }

  for(int i = _Y.size() - 2; i >= 0; i--) {
    if(ymax > _Y[i]) {
      _AOIMaxYi = i + 1;
      break;
    }
  }

  return retval;
}

PetscErrorCode CoarseGrid::Interpolate(const std::string& varname,
                                       double x, double y, double z, double t,
                                       double* value) {
  PetscErrorCode retval = 0;

  std::map<std::string, double*>::const_iterator mapiter;
  std::vector<std::string> dims;
  std::vector<unsigned int> start, count;
  double* buf = 0;
  size_t buflen = 0;


  retval = _Pio->inq_vardims(varname, dims);  CHKERRQ(retval);

  mapiter = _VarCache.find(varname);
  if(mapiter == _VarCache.end()) {
    // Then the variable has not yet been loaded and cached

    // TODO: this lazy-load approach might cause problems since I/O calls may
    // be collective...  For example, a process that doesn't have any values
    // in the NMS might never try to do an interpolation.

    buflen = 1;
    start.reserve(dims.size());
    count.reserve(dims.size());
    for(size_t i = 0; i < dims.size(); i++) {
      if(dims[i] == "x") {
        start[i] = _AOIMinXi;
        count[i] = _AOIMaxXi - _AOIMinXi + 1;
      } else if(dims[i] == "y") {
        start[i] = _AOIMinYi;
        count[i] = _AOIMaxYi - _AOIMinYi + 1;
      } else if(dims[i] == "z") {
        start[i] = 0;
        count[i] = _Z.size();
      } else if(dims[i] == "t") {
        start[i] = 0;
        count[i] = _T.size();;
      }

      buflen *= count[i];
    }

    buf = new double[buflen];
    _Pio->get_vara_double(varname, start, count, buf);
    _VarCache[varname] = buf;

  } else {
    buf = mapiter->second;
  }

  if(buf && dims.size() == 3) {
    // Assume t,x,y

    int xi1, xi2, yi1, yi2, ti1, ti2;
    int localw, localh;
    double x1, x2, y1, y2, t1, t2;
    double denom;
    double weights[8], values[8];

    for(int i = _AOIMinXi + 1; i <= _AOIMaxXi; i++) {
      if(_X[i - 1] <= x && x <= _X[i]) {
        xi1 = i - 1;
        xi2 = i;
        x1 = _X[xi1];
        x2 = _X[xi2];
        break;
      }
    }

    for(int i = _AOIMinYi + 1; i <= _AOIMaxYi; i++) {
      if(_Y[i - 1] <= y && y <= _Y[i]) {
        yi1 = i - 1;
        yi2 = i;
        y1 = _Y[yi1];
        y2 = _Y[yi2];
        break;
      }
    }

    for(size_t i = 0; i <= _T.size(); i++) {
      if(_T[i - 1] <= t && t <= _T[i]) {
        ti1 = i - 1;
        ti2 = i;
        t1 = _T[ti1];
        t2 = _T[ti2];
        break;
      }
    }

    denom = (t2 - t1)*(x2 - x1)*(y2 - y1);
    localw = _AOIMaxXi - _AOIMinXi + 1;
    localh = _AOIMaxYi - _AOIMinYi + 1;

    values[0] = buf[ti1*localw*localh + yi1*localw + xi1];
    weights[0] = (t2 - t)*(x2 - x)*(y2 - y) / denom;

    values[1] = buf[ti1*localw*localh + yi1*localw + xi2];
    weights[1] = (t2 - t)*(x - x1)*(y2 - y) / denom;

    values[2] = buf[ti1*localw*localh + yi2*localw + xi1];
    weights[2] = (t2 - t)*(x2 - x)*(y - y1) / denom;

    values[3] = buf[ti1*localw*localh + yi2*localw + xi2];
    weights[3] = (t2 - t)*(x - x1)*(y - y1) / denom;

    values[4] = buf[ti2*localw*localh + yi1*localw + xi1];
    weights[4] = (t - t1)*(x2 - x)*(y2 - y) / denom;

    values[5] = buf[ti2*localw*localh + yi1*localw + xi2];
    weights[5] = (t - t1)*(x - x1)*(y2 - y) / denom;

    values[6] = buf[ti2*localw*localh + yi2*localw + xi1];
    weights[6] = (t - t1)*(x2 - x)*(y - y1) / denom;

    values[7] = buf[ti2*localw*localh + yi2*localw + xi2];
    weights[7] = (t - t1)*(x - x1)*(y - y1) / denom;

    *value = 0.0;
    for(int i = 0; i < 8; i++) {
      *value += values[i] * weights[i];
    }

  } else if(buf && dims.size() == 4) {
    // Assume t,x,y,z
    fprintf(stderr, "4d interpolation not yet implemented.\n");
    retval = -1;
  }

  return retval;
}
