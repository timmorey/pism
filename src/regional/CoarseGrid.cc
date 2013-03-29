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
    _Filename(filename) {

  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  unsigned int xlen, ylen, zlen, tlen;
  PIO* pio = 0;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  pio = new PIO(comm, rank, "netcdf3");
  if(0 != pio->open(_Filename, PISM_NOWRITE)) {
    fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to open '%s'.\n", _Filename.c_str());
    delete pio;
    pio = 0;
  }

  if(pio) {
    if(0 != pio->inq_dimlen("x", xlen) ||
       0 != pio->inq_dimlen("y", ylen) ||
       0 != pio->inq_dimlen("z", zlen) ||
       0 != pio->inq_dimlen("time", tlen)) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to retrieve dimlen.\n");
      pio->close();
      delete pio;
      pio = 0;
    }
  }

  if(pio) {
    _X.reserve(xlen);
    _Y.reserve(ylen);
    _Z.reserve(zlen);
    _T.reserve(tlen);

    _AOIMaxXi = xlen - 1;
    _AOIMaxYi = ylen - 1;

    if((xlen > 0 && 0 != pio->get_1d_var("x", 0, xlen, _X)) ||
       (ylen > 0 && 0 != pio->get_1d_var("y", 0, ylen, _Y)) ||
       (zlen > 0 && 0 != pio->get_1d_var("z", 0, zlen, _Z)) ||
       (tlen > 0 && 0 != pio->get_1d_var("time", 0, tlen, _T))) {
      fprintf(stderr, "CoarseGrid::CoarseGrid: Failed to read dimension scales.\n");
      pio->close();
      delete pio;
      pio = 0;
    }
  }

  if(pio) {
    pio->close();
    delete pio;
    pio = 0;
  }
}

CoarseGrid::~CoarseGrid() {
  map<std::string, double*>::iterator mapiter;
  for(mapiter = _VarCache.begin(); mapiter != _VarCache.end(); mapiter++) {
    delete [] mapiter->second;
  }
}

PetscErrorCode CoarseGrid::SetAreaOfInterest(PetscReal xmin, PetscReal xmax,
                                             PetscReal ymin, PetscReal ymax) {
  PetscErrorCode retval = 0;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //printf("Rank %03d: CoarseGrid::SetAreaOfInterest(%f, %f, %f, %f)\n",
  //       rank, xmin, xmax, ymin, ymax);

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

  //printf("Rank %03d: CoarseGrid AOI: %d, %d, %d, %d\n",
  //      rank, _AOIMinXi, _AOIMaxXi, _AOIMinYi, _AOIMaxYi);

  return retval;
}

PetscErrorCode CoarseGrid::CacheVars(const std::vector<std::string>& varnames) {
  PetscErrorCode retval = 0;
  PIO* pio = 0;
  MPI_Comm comm = PETSC_COMM_WORLD;
  int rank;
  std::vector<std::string>::const_iterator nameiter;
  std::map<std::string, double*>::const_iterator mapiter;

  MPI_Comm_rank(comm, &rank);
  pio = new PIO(comm, rank, "netcdf3");
  if(0 != pio->open(_Filename, PISM_NOWRITE)) {
    fprintf(stderr, "CoarseGrid::CacheVars: Failed to open '%s'.\n", _Filename.c_str());
    delete pio;
    pio = 0;
  }

  if(pio) {
    for(nameiter = varnames.begin(); nameiter != varnames.end(); nameiter++) {
      mapiter = _VarCache.find(*nameiter);
      if(mapiter == _VarCache.end()) {
      // Then the variable has not yet been loaded and cached
        
        std::vector<std::string> dims;
        std::vector<unsigned int> start, count;
        double* buf = 0;
        size_t buflen = 0;
        
        retval = pio->inq_vardims(*nameiter, dims);  CHKERRQ(retval);
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
        pio->get_vara_double(*nameiter, start, count, buf);
        _VarCache[*nameiter] = buf;
        _VarDims[*nameiter] = dims;
      }
    }
  }

  if(pio) {
    pio->close();
    delete pio;
    pio = 0;
  }

  return retval;
}

PetscErrorCode CoarseGrid::Interpolate(const std::string& varname,
                                       double x, double y, double z, double t,
                                       double* value) {
  PetscErrorCode retval = 0;
//  int rank;

//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//  printf("Rank %03d: CoarseGrid::Interpolate('%s', %f, %f, %f, %f)\n",
//         rank, varname.c_str(), x, y, z, t);

  // NOTE: We don't want to do any _Pio-> calls in here, since those are
  // generally collective, and this function is not called equally at all 
  // processes.

  std::map<std::string, double*>::const_iterator mapiter;
  std::vector<std::string> dims;
  double* buf = 0;

  int ti1, ti2, xi1, xi2, yi1, yi2, zi1, zi2;
  int localw, localh, locald;
  double t1, t2, x1, x2, y1, y2, z1, z2;
  double denom;

  mapiter = _VarCache.find(varname);
  if(mapiter == _VarCache.end()) {
    fprintf(stderr, "Cannot interpolate '%s' - variable has not been loaded.\n", varname.c_str());
    retval = -1;
  } else {
    buf = mapiter->second;
    dims = _VarDims[varname];
  }

  if(buf) {
    if(dims.size() >= 1) {
      // Assume we have at least the t dim

      // Avoid seg fault in the case that we try to interpolate outside of the
      // available time range.  If everything is working right, we should only
      // ever get very slightly outside of the available time range.
      if(t < _T[0])
        t = _T[0];
      else if(t > _T[_T.size() - 1])
        t = _T[_T.size() - 1];

      for(size_t i = 1; i < _T.size(); i++) {
        if(_T[i - 1] <= t && t <= _T[i]) {
          ti1 = i - 1;
          ti2 = i;
          t1 = _T[ti1];
          t2 = _T[ti2];
          break;
        }
      }
    }
    
    if(buf && dims.size() >= 2) {
      // Assume we have at least the t,x dims
      
      localw = _AOIMaxXi - _AOIMinXi + 1;
      for(int i = _AOIMinXi + 1; i <= _AOIMaxXi; i++) {
        if(_X[i - 1] <= x && x <= _X[i]) {
          x1 = _X[i-1];
          x2 = _X[i];
          xi1 = i - 1 - _AOIMinXi;
          xi2 = i - _AOIMinXi;
          break;
        }
      }
    }
    
    if(buf && dims.size() >= 3) {
      // Assume we have at least the t,x,y dims
      
      localh = _AOIMaxYi - _AOIMinYi + 1;
      for(int i = _AOIMinYi + 1; i <= _AOIMaxYi; i++) {
        if(_Y[i - 1] <= y && y <= _Y[i]) {
          y1 = _Y[i-1];
          y2 = _Y[i];
          yi1 = i - 1 - _AOIMinYi;
          yi2 = i - _AOIMinYi;
          break;
        }
      }
    }
    
    if(buf && dims.size() >= 4) {
      // Assume we have at least the t,x,y,z dims
      
      locald = _Z.size();
      for(size_t i = 1; i < _Z.size(); i++) {
        if(_Z[i - 1] <= z && z <= _Z[i]) {
          zi1 = i - 1;
          zi2 = i;
          z1 = _Z[zi1];
          z2 = _Z[zi2];
          break;
        }
      }
    }
    
    if(buf && dims.size() == 3) {
      // Assume t,x,y
      int npts = 8;
      double weights[8], values[8];
      
      denom = (t2 - t1)*(x2 - x1)*(y2 - y1);
      
      values[0] = buf[ti1*localw*localh + xi1*localh + yi1];
      weights[0] = (t2 - t)*(x2 - x)*(y2 - y) / denom;
      
      values[1] = buf[ti1*localw*localh + xi2*localh + yi1];
      weights[1] = (t2 - t)*(x - x1)*(y2 - y) / denom;

      values[2] = buf[ti1*localw*localh + xi1*localh + yi2];
      weights[2] = (t2 - t)*(x2 - x)*(y - y1) / denom;
      
      values[3] = buf[ti1*localw*localh + xi2*localh + yi2];
      weights[3] = (t2 - t)*(x - x1)*(y - y1) / denom;
      
      values[4] = buf[ti2*localw*localh + xi1*localh + yi1];
      weights[4] = (t - t1)*(x2 - x)*(y2 - y) / denom;
      
      values[5] = buf[ti2*localw*localh + xi2*localh + yi1];
      weights[5] = (t - t1)*(x - x1)*(y2 - y) / denom;

      values[6] = buf[ti2*localw*localh + xi1*localh + yi2];
      weights[6] = (t - t1)*(x2 - x)*(y - y1) / denom;
      
      values[7] = buf[ti2*localw*localh + xi2*localh + yi2];
      weights[7] = (t - t1)*(x - x1)*(y - y1) / denom;

      *value = 0.0;
      for(int i = 0; i < npts; i++) {
        *value += values[i] * weights[i];
      }

    } else if(buf && dims.size() == 4) {
      // Assume t,x,y,z
      int npts = 16;
      double weights[16], values[16];
      int tstride = localw*localh*locald;
      int xstride = localh*locald;
      int ystride = locald;
      
      denom = (t2 - t1)*(x2 - x1)*(y2 - y1)*(z2 - z1);

      values[0] = buf[ti1*tstride + xi1*xstride + yi1*ystride + zi1];
      weights[0] = (t2 - t)*(x2 - x)*(y2 - y)*(z2 - z) / denom;

      values[1] = buf[ti1*tstride + xi1*xstride + yi1*ystride + zi2];
      weights[1] = (t2 - t)*(x2 - x)*(y2 - y)*(z - z1) / denom;

      values[2] = buf[ti1*tstride + xi1*xstride + yi2*ystride + zi1];
      weights[2] = (t2 - t)*(x2 - x)*(y - y1)*(z2 - z) / denom;

      values[3] = buf[ti1*tstride + xi1*xstride + yi2*ystride + zi2];
      weights[3] = (t2 - t)*(x2 - x)*(y - y1)*(z - z1) / denom;

      values[4] = buf[ti1*tstride + xi2*xstride + yi1*ystride + zi1];
      weights[4] = (t2 - t)*(x - x1)*(y2 - y)*(z2 - z) / denom;

      values[5] = buf[ti1*tstride + xi2*xstride + yi1*ystride + zi2];
      weights[5] = (t2 - t)*(x - x1)*(y2 - y)*(z - z1) / denom;

      values[6] = buf[ti1*tstride + xi2*xstride + yi2*ystride + zi1];
      weights[6] = (t2 - t)*(x - x1)*(y - y1)*(z2 - z) / denom;

      values[7] = buf[ti1*tstride + xi2*xstride + yi2*ystride + zi2];
      weights[7] = (t2 - t)*(x - x1)*(y - y1)*(z - z1) / denom;

      values[8] = buf[ti2*tstride + xi1*xstride + yi1*ystride + zi1];
      weights[8] = (t - t1)*(x2 - x)*(y2 - y)*(z2 - z) / denom;

      values[9] = buf[ti2*tstride + xi1*xstride + yi1*ystride + zi2];
      weights[9] = (t - t1)*(x2 - x)*(y2 - y)*(z - z1) / denom;

      values[10] = buf[ti2*tstride + xi1*xstride + yi2*ystride + zi1];
      weights[10] = (t - t1)*(x2 - x)*(y - y1)*(z2 - z) / denom;

      values[11] = buf[ti2*tstride + xi1*xstride + yi2*ystride + zi2];
      weights[11] = (t - t1)*(x2 - x)*(y - y1)*(z - z1) / denom;

      values[12] = buf[ti2*tstride + xi2*xstride + yi1*ystride + zi1];
      weights[12] = (t - t1)*(x - x1)*(y2 - y)*(z2 - z) / denom;

      values[13] = buf[ti2*tstride + xi2*xstride + yi1*ystride + zi2];
      weights[13] = (t - t1)*(x - x1)*(y2 - y)*(z - z1) / denom;

      values[14] = buf[ti2*tstride + xi2*xstride + yi2*ystride + zi1];
      weights[14] = (t - t1)*(x - x1)*(y - y1)*(z2 - z) / denom;

      values[15] = buf[ti2*tstride + xi2*xstride + yi2*ystride + zi2];
      weights[15] = (t - t1)*(x - x1)*(y - y1)*(z - z1) / denom;

      *value = 0.0;
      for(int i = 0; i < npts; i++) {
        *value += values[i] * weights[i];
      }

    } else {
      fprintf(stderr, "%dd interpolation not yet implemented.\n", (int)dims.size());
      retval = -1;
    }
  }
  
  return retval;
}
