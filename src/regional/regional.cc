// Copyright (C) 2010, 2011, 2012 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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

#include "CoarseGrid.hh"
#include "regional.hh"

SIAFD_Regional::SIAFD_Regional(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &c,
                               CoarseGrid* cg)
  : SIAFD(g, e, c),
    coarse_grid(cg) {

}

PetscErrorCode SIAFD_Regional::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SIAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SIA solver...\n"); CHKERRQ(ierr);

  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(grid.com, 1, "no_model_mask is not available");

  usurfstore = dynamic_cast<IceModelVec2S*>(vars.get("usurfstore"));
  if (usurfstore == NULL) SETERRQ(grid.com, 1, "usurfstore is not available");

  return 0;
}

PetscErrorCode SIAFD_Regional::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  ierr = SIAFD::compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);

  IceModelVec2Int nmm = *no_model_mask;
  IceModelVec2S hst = *usurfstore; // convenience

  const int Mx = grid.Mx, My = grid.My;
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = nmm.begin_access(); CHKERRQ(ierr);
  ierr = hst.begin_access(); CHKERRQ(ierr);

  // TODO: we are ignoring the zero_grad_where_no_model flag here, if we have a
  // CoarseGrid to interpolate NMS values.

  double t = grid.time->current();
  PetscInt GHOSTS = 1;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {

      // x-component, i-offset
      if (nmm(i, j) > 0.5 || nmm(i + 1, j) > 0.5) {
        if (i < 0 || i + 1 > Mx - 1) {
          h_x(i, j, 0) = 0.0;
        } else {
          if(coarse_grid) {
            double hst_ip1_j, hst_i_j;
            coarse_grid->Interpolate("usurf", grid.x[i+1], grid.y[j], 0.0, t, &hst_ip1_j);
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j], 0.0, t, &hst_i_j);
            //printf("Interpolated usurf: %f, %f\n", hst_ip1_j, hst_i_j);
            h_x(i, j, 0) = (hst_ip1_j - hst_i_j) / dx;            
          } else {
            h_x(i, j, 0) = (hst(i + 1, j) - hst(i, j)) / dx;
          }
        }
      }

      // x-component, j-offset
      if (nmm(i - 1, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
          nmm(i - 1, j)     > 0.5 || nmm(i + 1, j)     > 0.5) {
        if (i - 1 < 0 || j + 1 > My - 1 || i + 1 > Mx - 1) {
          h_x(i, j, 1) = 0.0;
        } else {
          if(coarse_grid) {
            double hst_ip1_jp1, hst_ip1_j, hst_im1_jp1, hst_im1_j;
            coarse_grid->Interpolate("usurf", grid.x[i+1], grid.y[j+1], 0.0, t, &hst_ip1_jp1);
            coarse_grid->Interpolate("usurf", grid.x[i+1], grid.y[j], 0.0, t, &hst_ip1_j);
            coarse_grid->Interpolate("usurf", grid.x[i-1], grid.y[j+1], 0.0, t, &hst_im1_jp1);
            coarse_grid->Interpolate("usurf", grid.x[i-1], grid.y[j], 0.0, t, &hst_im1_j);
            //printf("Interpolated usurf: %f, %f, %f, %f\n", hst_ip1_jp1, hst_ip1_j, hst_im1_jp1, hst_im1_j);
            h_x(i, j, 1) = ( + hst_ip1_jp1 + hst_ip1_j
                             - hst_im1_jp1 - hst_im1_j ) / (4.0 * dx);
          } else {
            h_x(i, j, 1) = ( + hst(i + 1, j + 1) + hst(i + 1, j)
                             - hst(i - 1, j + 1) - hst(i - 1, j) ) / (4.0 * dx);
          }
        }
      }

      // y-component, i-offset
      if (nmm(i, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
          nmm(i, j - 1) > 0.5 || nmm(i + 1, j - 1) > 0.5) {
        if (i < 0 || j + 1 > My - 1 || i + 1 > Mx - 1 || j - 1 < 0) {
          h_y(i, j, 0) = 0.0;
        } else {
          if(coarse_grid) {
            double hst_ip1_jp1, hst_i_jp1, hst_ip1_jm1, hst_i_jm1;
            coarse_grid->Interpolate("usurf", grid.x[i+1], grid.y[j+1], 0.0, t, &hst_ip1_jp1);
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j+1], 0.0, t, &hst_i_jp1);
            coarse_grid->Interpolate("usurf", grid.x[i+1], grid.y[j-1], 0.0, t, &hst_ip1_jm1);
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j-1], 0.0, t, &hst_i_jm1);
            //printf("Interpolated usurf: %f, %f, %f, %f\n", hst_ip1_jp1, hst_i_jp1, hst_ip1_jm1, hst_i_jm1);
            h_y(i, j, 0) = ( + hst_ip1_jp1 + hst_i_jp1
                             - hst_ip1_jm1 - hst_i_jm1 ) / (4.0 * dy);
          } else {
            h_y(i, j, 0) = ( + hst(i + 1, j + 1) + hst(i, j + 1)
                             - hst(i + 1, j - 1) - hst(i, j - 1) ) / (4.0 * dy);
          }
        }
      }

      // y-component, j-offset
      if (nmm(i, j) > 0.5 || nmm(i, j + 1) > 0.5) {
        if (j < 0 || j + 1 > My - 1) {
          h_y(i, j, 1) = 0.0;
        } else {
          if(coarse_grid) {
            double hst_i_jp1, hst_i_j;
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j+1], 0.0, t, &hst_i_jp1);
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j], 0.0, t, &hst_i_j);
            //printf("Interpolated usurf: %f, %f\n", hst_i_jp1, hst_i_j);
            h_y(i, j, 1) = (hst_i_jp1 - hst_i_j) / dy;
          } else {
            h_y(i, j, 1) = (hst(i, j + 1) - hst(i, j)) / dy;
          }
        }
      }
    }
  }
  ierr = nmm.end_access(); CHKERRQ(ierr);
  ierr = hst.end_access(); CHKERRQ(ierr);
  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}

SSAFD_Regional::SSAFD_Regional(IceGrid &g, IceBasalResistancePlasticLaw &b, EnthalpyConverter &e,
                               const NCConfigVariable &c, CoarseGrid* cg)
  : SSAFD(g, b, e, c),
    coarse_grid(cg) {

}

PetscErrorCode SSAFD_Regional::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SSA solver...\n"); CHKERRQ(ierr);

  if (config.get_flag("ssa_dirichlet_bc")) {
    ierr = verbPrintf(2,grid.com,"  using stored SSA velocities as Dirichlet B.C. in the no_model_strip...\n"); 
    CHKERRQ(ierr);
  }

  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(grid.com, 1, "no_model_mask is not available");

  usurfstore = dynamic_cast<IceModelVec2S*>(vars.get("usurfstore"));
  if (usurfstore == NULL) SETERRQ(grid.com, 1, "usurfstore is not available");

  thkstore = dynamic_cast<IceModelVec2S*>(vars.get("thkstore"));
  if (thkstore == NULL) SETERRQ(grid.com, 1, "thkstore is not available");

  return 0;
}

PetscErrorCode SSAFD_Regional::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  ierr = SSAFD::compute_driving_stress(result); CHKERRQ(ierr);

  const PetscReal standard_gravity = config.get("standard_gravity"),
    ice_rho = config.get("ice_density");
  IceModelVec2Int nmm = *no_model_mask;

  // TODO: we are ignoring the zero_grad_where_no_model flag here if we have a
  // CoarseGrid to interpolate NMS values.

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = nmm.begin_access(); CHKERRQ(ierr);
  ierr = usurfstore->begin_access(); CHKERRQ(ierr);
  ierr = thkstore->begin_access(); CHKERRQ(ierr);
  double t = grid.time->current();
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar pressure = ice_rho * standard_gravity * (*thkstore)(i,j);
      if (pressure <= 0) pressure = 0;

      if (nmm(i, j) > 0.5 || nmm(i - 1, j) > 0.5 || nmm(i + 1, j) > 0.5) {
        if (i - 1 < 0 || i + 1 > grid.Mx - 1) {
          result(i, j).u = 0;
        } else {
          if(coarse_grid) {
            double x1, x2, dx;
            coarse_grid->Interpolate("usurf", grid.x[i-1], grid.y[j], 0.0, t, &x1);
            coarse_grid->Interpolate("usurf", grid.x[i+1], grid.y[j], 0.0, t, &x2);
            dx = grid.dx;
            //printf("Interpolated usurf: %f, %f\n", x1, x2);
            //printf("Interpolated diff_x: %f\n", (x2 - x1) / (2.0 * dx));
            result(i, j).u = - pressure * ((x2 - x1) / (2.0 * dx));
          } else {
            result(i, j).u = - pressure * usurfstore->diff_x(i,j);
          }
        }
      }

      if (nmm(i, j) > 0.5 || nmm(i, j - 1) > 0.5 || nmm(i, j + 1) > 0.5) {
        if (j - 1 < 0 || j + 1 > grid.My - 1) {
          result(i, j).v = 0;
        } else {
          if(coarse_grid) {
            double y1, y2, dy;
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j-1], 0.0, t, &y1);
            coarse_grid->Interpolate("usurf", grid.x[i], grid.y[j+1], 0.0, t, &y2);
            dy = grid.dy;
            //printf("Interpolated usurf: %f, %f\n", y1, y2);
            //printf("Interpolated diff_y: %f\n", (y2 - y1) / (2.0 * dy));
            result(i, j).v = - pressure * ((y2 - y1) / (2.0 * dy));
          } else {
            result(i, j).v = - pressure * usurfstore->diff_y(i,j);
          }
        }
      }

    }
  }
  ierr = usurfstore->end_access(); CHKERRQ(ierr);
  ierr = thkstore->end_access(); CHKERRQ(ierr);
  ierr = nmm.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMRegionalDefaultYieldStress::init(PISMVars &vars) {
  PetscErrorCode ierr;
  PetscInt v = getVerbosityLevel(); // turn off second, redundant init message
  ierr = setVerbosityLevel(1); CHKERRQ(ierr);
  ierr = PISMMohrCoulombYieldStress::init(vars); CHKERRQ(ierr);
  ierr = setVerbosityLevel(v); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
    "  using the regional version with strong till in no_model_mask==1 area ...\n");
    CHKERRQ(ierr);
  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(grid.com, 1, "no_model_mask is not available");
  return 0;
}


PetscErrorCode PISMRegionalDefaultYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  PetscErrorCode ierr;
  
  // do whatever you normally do
  ierr = PISMMohrCoulombYieldStress::basal_material_yield_stress(result); CHKERRQ(ierr);

  // now set result=tauc to a big value in no_model_strip
  ierr = no_model_mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if ((*no_model_mask)(i,j) > 0.5) {
        result(i,j) = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar
      }
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = no_model_mask->end_access(); CHKERRQ(ierr);
  return 0;
}

