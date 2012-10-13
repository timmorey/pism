// Copyright (C) 2011 Matthias Mengel and Torsten Albrecht
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

#include <cmath>
#include <cstring>
#include <petscdmda.h>

#include "iceModel.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "PISMOcean.hh"
#include "pism_options.hh"


//! \file iMpartgridground.cc Methods implementing PIK option -part_grid_ground
//! matthias.mengel@pik

PetscReal IceModel::get_average_thickness_fg(planeStar<int> M, planeStar<PetscScalar> H, planeStar<PetscScalar> h,
                                             planeStar<PetscScalar> Q, planeStar<PetscScalar> Qssa, PetscReal bed_ij,
                                             PetscScalar &sia_ssa_coeff) {
  PetscErrorCode ierr;
  bool margin_coeff_set;
  PetscReal margin_coeff;
  ierr = PISMOptionsReal("-margin_coeff_ground", "specifies the ratio of partially filled grid box height to its surrounding neighbours", margin_coeff,  margin_coeff_set); CHKERRQ(ierr);
  if (!margin_coeff_set) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Please specify scale coefficient for scaleMargin method.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  const bool const_coeff = config.get_flag("do_const_pgg_scaling");

  Mask m;
  PetscInt N = 0;
  PetscReal H_average = 0.0;

  // calculate Href as average height of the neighbouring grounded grid cells
  // get mean ice thickness over adjacent grounded ice shelf neighbors

  // if a direct neighbour is grounded ice, then add H to H_average and increase N.
  if (m.grounded_ice(M.e)) { H_average += H.e; N++; }
  if (m.grounded_ice(M.w)) { H_average += H.w; N++; }
  if (m.grounded_ice(M.n)) { H_average += H.n; N++; }
  if (m.grounded_ice(M.s)) { H_average += H.s; N++; }

  if (N > 0) {
    // average by division by grounded neighbours
    H_average = H_average / N;
  } else {
    // assume constant H_average for isolated cells
    H_average = 200.0;
  }



  if( const_coeff ){
    PetscSynchronizedPrintf(grid.com,"const pgg coeff=%f, Havg=%e\n",margin_coeff, H_average);
    return H_average * margin_coeff;
  }
  //PetscSynchronizedPrintf(grid.com,"not printed if pgg simple on\n");

  //
  //if( m.grounded_ice(M.ij) ){
  PetscReal h_average = 0.0;
  PetscReal H_average_FromBed = 0.0;
  N = 0;
  if (m.grounded_ice(M.e)) { h_average += h.e; N++; }
  if (m.grounded_ice(M.w)) { h_average += h.w; N++; }
  if (m.grounded_ice(M.n)) { h_average += h.n; N++; }
  if (m.grounded_ice(M.s)) { h_average += h.s; N++; }

  // ice thickness is difference between elevation and bed
  if (N > 0) {
    H_average_FromBed = h_average / N - bed_ij;
  } else {
    // assume constant H_average for isolated cells
    H_average_FromBed = 200.0 - bed_ij;
  }

  // decide on which Href to use, this is usually H_average for downward sloping and H_average_FromBed for upward sloping ground.
  PetscSynchronizedPrintf(grid.com,"Havg=%e, Havg_bed=%e, bed=%e\n",H_average, H_average_FromBed,bed_ij);
  H_average = PetscMin(H_average, H_average_FromBed);
  //}

  PetscReal Qssa_max      = 0.0;
  PetscReal sum_max       = 0.0;

  if( sum_max < PetscAbs(Q.e + Qssa.e) ){ sum_max = PetscAbs(Q.e + Qssa.e); Qssa_max = PetscAbs(Qssa.e);}
  if( sum_max < PetscAbs(Q.w + Qssa.w) ){ sum_max = PetscAbs(Q.w + Qssa.w); Qssa_max = PetscAbs(Qssa.w);}
  if( sum_max < PetscAbs(Q.n + Qssa.n) ){ sum_max = PetscAbs(Q.n + Qssa.n); Qssa_max = PetscAbs(Qssa.n);}
  if( sum_max < PetscAbs(Q.s + Qssa.s) ){ sum_max = PetscAbs(Q.s + Qssa.s); Qssa_max = PetscAbs(Qssa.s);}

  if (sum_max > 0.0) {
    sia_ssa_coeff = margin_coeff + PetscMin(Qssa_max/sum_max,1.0)*(1.0-margin_coeff);
  } else {
    sia_ssa_coeff = 1.0;
  }

  return H_average * sia_ssa_coeff;
}

PetscErrorCode IceModel::killLonelyPGGCells() {
  PetscErrorCode ierr;
  // looking for one-grid-cell partially filled grid cells, that have 4 neighbors of thickness H=0
  const bool vpik = config.get_flag("verbose_pik_messages");

  PetscReal sea_level;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com,2, "PISM ERROR: ocean == NULL");
  }
  double ocean_rho = config.get("sea_water_density"),
         ice_rho   = config.get("ice_density");

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vHrefGround.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);

  PetscReal C = (1.0 - ice_rho / ocean_rho);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      planeStar<PetscScalar> thk = vH.star(i, j),
                             bed = vbed.star(i, j);

      bool all_4neighbors_ungrounded =
        (bed.e + thk.e <  sea_level + C * thk.e) &&
        (bed.w + thk.w <  sea_level + C * thk.w) &&
        (bed.n + thk.n <  sea_level + C * thk.n) &&
        (bed.s + thk.s <  sea_level + C * thk.s);

      if ( vHrefGround(i, j) > 0.0 && all_4neighbors_ungrounded) {
        vHrefGround(i, j) = 0.0;
        if (vpik) {
              PetscSynchronizedPrintf(grid.com,
                "PISM-PIK INFO: [rank %d] killed lonely PGG cell at i = %d, j = %d\n",
                grid.rank, i, j);
        }
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vHrefGround.end_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  return 0;
}
