// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

static char help[] =
  "Ice sheet driver for PISM ice sheet simulations, initialized from data.\n"
  "The basic PISM executable for evolution runs.\n";

#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"

#include "pism_options.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"
#include "PISMProf.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  PISMLogStagePush(PISM_INIT_STAGE);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMR %s (basic evolution run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    ierr = check_old_option_and_stop(com, "-boot_from", "-boot_file"); CHKERRQ(ierr); 

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
    string usage =
      "  pismr {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismr", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pismr", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides, true); CHKERRQ(ierr);

    IceGrid g(com, rank, size, config);
    IceModel m(g, config, overrides);

    // Initialize boundary models:
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;

    POFactory po(g, config);
    PISMOceanModel *ocean;

    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    surface->attach_atmosphere_model(atmosphere);

    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);
    ierr = m.setExecName("pismr"); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    PISMLogStagePop(); // PISM_INIT_STAGE

    PISMLogStagePush(PISM_STEPPING_STAGE);

    ierr = m.run(); CHKERRQ(ierr);

    PISMLogStagePop(); // PISM_STEPPING_STAGE

    PISMLogStagePush(PISM_PRIMARY_OUTPUT_STAGE);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);
    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamed.nc"); CHKERRQ(ierr);

    PISMLogStagePop(); // PISM_PRIMARY_OUTPUT_STAGE

    PISMLogStagePush(PISM_SHUTDOWN_STAGE);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
