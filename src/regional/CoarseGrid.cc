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

  _Pio = new PIO(comm, rank, "netcdf3");
  if(0 != _Pio->open(filename, PISM_NOWRITE)) {
    fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to open '%s'.\n", filename.c_str());
    delete _Pio;
    _Pio = 0;
  }

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

  if(_Pio) {
    _X.reserve(xlen);
    _Y.reserve(ylen);
    _Z.reserve(zlen);
    _T.reserve(tlen);

    _AOIMaxXi = xlen - 1;
    _AOIMaxYi = ylen - 1;

    if((xlen > 0 && 0 != _Pio->get_1d_var("x", 0, xlen, _X)) ||
       (ylen > 0 && 0 != _Pio->get_1d_var("y", 0, ylen, _Y)) ||
       (zlen > 0 && 0 != _Pio->get_1d_var("z", 0, zlen, _Z)) ||
       (tlen > 0 && 0 != _Pio->get_1d_var("time", 0, tlen, _T))) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to read dimension scales.\n");
      _Pio->close();
      delete _Pio;
      _Pio = 0;
    }
  }
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

PetscErrorCode CoarseGrid::CacheVars(const std::vector<std::string>& varnames) {
  PetscErrorCode retval = 0;

  std::vector<std::string>::const_iterator nameiter;
  std::map<std::string, double*>::const_iterator mapiter;
  std::vector<std::string> dims;
  std::vector<unsigned int> start, count;
  double* buf = 0;
  size_t buflen = 0;

  for(nameiter = varnames.begin(); nameiter != varnames.end(); nameiter++) {
    mapiter = _VarCache.find(*nameiter);
    if(mapiter == _VarCache.end()) {
      // Then the variable has not yet been loaded and cached

      retval = _Pio->inq_vardims(*nameiter, dims);  CHKERRQ(retval);
      buflen = 1;
      for(size_t i = 0; i < dims.size(); i++) {
        if(dims[i] == "x") {
          start.push_back(_AOIMinXi);
          count.push_back(_AOIMaxXi - _AOIMinXi + 1);
        } else if(dims[i] == "y") {
          start.push_back(_AOIMinYi);
          count.push_back(_AOIMaxYi - _AOIMinYi + 1);
        } else if(dims[i] == "z") {
          start.push_back(0);
          count.push_back(_Z.size());
        } else if(dims[i] == "time") {
          start.push_back(0);
          count.push_back(_T.size());
        }
        
        buflen *= count[i];
      }

      buf = new double[buflen];
      _Pio->get_vara_double(*nameiter, start, count, buf);
      _VarCache[*nameiter] = buf;
      _VarDims[*nameiter] = dims;
    }
  }

  return retval;
}

PetscErrorCode CoarseGrid::Interpolate(const std::string& varname,
                                       double x, double y, double z, double t,
                                       double* value) {
  PetscErrorCode retval = 0;

  // NOTE: We don't want to do any _Pio-> calls in here, since those are
  // generally collective, and this function is not called equally at all 
  // processes.

  std::map<std::string, double*>::const_iterator mapiter;
  std::vector<std::string> dims;
  double* buf = 0;

  mapiter = _VarCache.find(varname);
  if(mapiter == _VarCache.end()) {
    fprintf(stderr, "Cannot interpolate '%s' - variable has not been loaded.\n", varname.c_str());
    retval = -1;
  } else {
    buf = mapiter->second;
    dims = _VarDims[varname];
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
        x1 = _X[i-1];
        x2 = _X[i];
        xi1 = i - 1 - _AOIMinXi;
        xi2 = i - _AOIMinXi;
        break;
      }
    }

    for(int i = _AOIMinYi + 1; i <= _AOIMaxYi; i++) {
      if(_Y[i - 1] <= y && y <= _Y[i]) {
        y1 = _Y[i-1];
        y2 = _Y[i];
        yi1 = i - 1 - _AOIMinYi;
        yi2 = i - _AOIMinYi;
        break;
      }
    }

    for(size_t i = 1; i < _T.size(); i++) {
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
