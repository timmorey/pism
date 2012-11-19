// Copyright (C) 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "PISMNCFile.hh"

#include <cstdio>               // fprintf, stderr

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

#include <boost/algorithm/string.hpp>

PISMNCFile::PISMNCFile(MPI_Comm c, int r)
  : rank(r), com(c), ncid(-1), define_mode(false) {

  MPI_Info_create(&mpi_info);
  init_hints();
}

PISMNCFile::~PISMNCFile() {
  MPI_Info_free(&mpi_info);
}

string PISMNCFile::get_filename() const {
  return filename;
}

int PISMNCFile::put_att_double(string variable_name, string att_name, PISM_IO_Type nctype, double value) const {
  vector<double> tmp(1);
  tmp[0] = value;
  return put_att_double(variable_name, att_name, nctype, tmp);
}

//! \brief Prints an error message; for debugging.
void PISMNCFile::check(int return_code) const {
  if (return_code != NC_NOERR) {
    fprintf(stderr, "NC_ERR: %s\n", nc_strerror(return_code));
  }
}

void PISMNCFile::init_hints() const {
  string env = "";
  string hint = "";
  char name[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
  char* temp = 0;
  size_t pos = 0;

  temp = getenv("PISM_MPIIO_HINTS");
  if(temp)
    env = temp;

  while(env.size() > 0) {
    pos = env.find(";");
    if(string::npos == pos) {
      hint = env;
      env = "";
    } else {
      hint = env.substr(0, pos);
      env = env.substr(pos + 1);
    }

    pos = hint.find("=");
    if(string::npos == pos) {
      if(0 == rank)
        fprintf(stderr, "Incorrect syntax in hint: '%s'\n", hint.c_str());
    } else {
      strcpy(name, hint.substr(0, pos).c_str());
      strcpy(value, hint.substr(pos + 1).c_str());
      if(0 == rank)
        printf("Setting hint %s=%s\n", name, value);
      MPI_Info_set(mpi_info, name, value);
    }
  }
}
