// Copyright (C) 2010, 2011, 2012, 2013 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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
  "Ice sheet driver for PISM regional (outlet glacier) simulations, initialized\n"
  "from data.\n";

#include <petsc.h>

#include "IceGrid.hh"
#include "IceRegionalModel.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"
#include "regional.hh"

//! \file pismo.cc A regional (outlet glacier) model form of PISM.
/*! \file pismo.cc 
The classes in this file modify basic PISM whole ice sheet modeling assumptions.
Normally in PISM the ice sheet occupies a continent which is surrounded by
ocean.  Or at least PISM assumes that the edge of the computational domain is in
a region with strong ablation that the ice will not cross.

Here, by contrast, we add a strip around the edge of the computational domain
(variable \c no_model_mask and option \c -no_model_strip).  Various
simplifications and boundary conditions are enforced in this script:
* the surface gradient computation is made trivial,
* the driving stress does not change during the run but instead comes from
the gradient of a saved surface elevation, and
* the base is made strong so that no sliding occurs.

Also options \c -force_to_thk and variable \c ftt_mask play a role in isolating
the modeled outlet glacier.  But there is no code here for that purpose. 
Instead see the PSForceThickness surface model modifier class.
 */

int main(int argc, char *argv[]) {
  char temp[256];
  PetscBool coarseGridFileSet;
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  MPI_Comm    com = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMO %s (regional outlet-glacier run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
    string usage =
      "  pismo {-i IN.nc|-boot_file IN.nc} [-no_model_strip X] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "  -no_model_strip X (re-)set width of no-model strip along edge of\n"
      "              computational domain to X km\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((!iset) && (!bfset)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismo", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;
      required.clear();
      ierr = show_usage_check_req_opts(com, "pismo", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides, true); CHKERRQ(ierr);

    ierr = PetscOptionsBegin(com, "", "Embedded Modeling Options", ""); CHKERRQ(ierr);
    ierr = PetscOptionsString("-coarse_grid_file",
                              "A file that contains a coarse grid suitable for "
                              "interpolating boundary conditions in the "
                              "regional model.",
                              "", "", temp, 256, &coarseGridFileSet);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // initialize the ice dynamics model
    IceGrid g(com, rank, size, config);
    IceRegionalModel m(g, config, overrides);
    ierr = m.setExecName("pismo"); CHKERRQ(ierr);

    // initialize boundary models
    // factories allow runtime choice
    PAFactory pa(g, config);
    PSFactory ps(g, config);
    POFactory po(g, config);
    // now read options and choose
    PISMAtmosphereModel *atmosphere;
    PISMSurfaceModel    *surface;
    PISMOceanModel      *ocean;
    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    surface->attach_atmosphere_model(atmosphere); // IceModel m does not see atmosphere
    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);

    ierr = m.init(); CHKERRQ(ierr);

    if(coarseGridFileSet) {
      printf("Loading coarse grid file '%s'...\n", temp);
      m.attach_coarse_grid(temp);
    }

    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamed_regional.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

