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
#include <petscda.h>

#include "iceModel.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "PISMOcean.hh"
// #include "pism_options.hh"


//! \file iMpartgridground.cc Methods implementing PIK option -part_grid_ground
//! matthias.mengel@pik

PetscReal IceModel::get_average_thickness_fg(planeStar<int> M, planeStar<PetscScalar> H, planeStar<PetscScalar> h,
                                             planeStar<PetscScalar> Q, planeStar<PetscScalar> Qssa, PetscReal bed_ij, PetscScalar &sia_ssa_coeff) {
  PetscErrorCode ierr;
  bool margin_coeff_set;
  PetscReal margin_coeff;
  ierr = PISMOptionsReal("-margin_coeff_ground", "specifies the ratio of partially filled grid box height to its surrounding neighbours", margin_coeff,  margin_coeff_set); CHKERRQ(ierr);
  if (!margin_coeff_set) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Please specify scale coefficient for scaleMargin method.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

//   ierr = vTestVar.begin_access(); CHKERRQ(ierr);

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

  if( m.grounded_ice(M.ij)){
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
    H_average = PetscMin(H_average, H_average_FromBed);
  }

  PetscReal Qssa_max      = 0.0;
  PetscReal sum_max       = 0.0;
//   PetscReal sia_ssa_coeff = 0.0;
  if( sum_max < PetscAbs(Q.e + Qssa.e) ) sum_max = PetscAbs(Q.e + Qssa.e); Qssa_max = PetscAbs(Qssa.e);
  if( sum_max < PetscAbs(Q.w + Qssa.w) ) sum_max = PetscAbs(Q.w + Qssa.w); Qssa_max = PetscAbs(Qssa.w);
  if( sum_max < PetscAbs(Q.n + Qssa.n) ) sum_max = PetscAbs(Q.n + Qssa.n); Qssa_max = PetscAbs(Qssa.n);
  if( sum_max < PetscAbs(Q.s + Qssa.s) ) sum_max = PetscAbs(Q.s + Qssa.s); Qssa_max = PetscAbs(Qssa.s);

  if (sum_max > 0) {
    sia_ssa_coeff = margin_coeff + PetscMin(Qssa_max/sum_max,1)*(1-margin_coeff);
  } else {
    sia_ssa_coeff = 1.0;
  }

//   vTestVar(i,j) = sia_ssa_coeff;
//   ierr = vTestVar.end_access(); CHKERRQ(ierr);
  return H_average * sia_ssa_coeff;
}

PetscErrorCode IceModel::killLonelyPGGCells() {
  PetscErrorCode ierr;
  // looking for one-grid-cell partially filled grid cells, that have 4 neighbors of thickness H=0
  const bool vpik = config.get_flag("verbose_pik_messages");

  PetscReal sea_level;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(2, "PISM ERROR: ocean == NULL"); }
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
//         vMask(i, j) = MASK_ICE_FREE_OCEAN;
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

PetscReal IceModel::get_average_thickness_g(planeStar<int> M, planeStar<PetscScalar> H, PetscReal bed_ij, PetscInt i, PetscInt j) {
  
  //////////// get part_grid_ground choice 
  PetscErrorCode ierr;
  bool calcMethod_set;
  string HrefCalcMethod;
  set<string> HrefCalc_choices;
  HrefCalc_choices.insert("Hav");
  HrefCalc_choices.insert("hav");
  HrefCalc_choices.insert("fixed");
  HrefCalc_choices.insert("scaleMargin");
  ierr = PISMOptionsList(grid.com,"-part_grid_ground", "specifies the Href calculation method",
         HrefCalc_choices, "Hav", HrefCalcMethod, calcMethod_set); CHKERRQ(ierr);
  if (!calcMethod_set) {
    ierr = PetscPrintf(grid.com,
           "PISM ERROR: Please specify an method to calc HrefGround.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }
  ////////////
  
  PetscInt N = 0;
  Mask m;
  PetscReal H_average = 0.0;
  PetscReal margin_coeff;
  ierr = vTestVar.begin_access(); CHKERRQ(ierr);

  if (HrefCalcMethod == "fixed"){

    bool fixedHeight_set;
    PetscReal fixedHeight;
    ierr = PISMOptionsReal("-fixedHeight", "specifies fixed Height for H_average",
                          fixedHeight,  fixedHeight_set); CHKERRQ(ierr);
    if (!fixedHeight_set) {
      ierr = PetscPrintf(grid.com,
            "PISM ERROR: Please specify a fixed height for H_average.\n");
      CHKERRQ(ierr);
      PISMEnd();
    }
    H_average = fixedHeight;
    return H_average;
  }

  if (HrefCalcMethod == "scaleMargin"){
    
    bool margin_coeff_set;
    ierr = PISMOptionsReal("-margin_coeff_ground", "specifies the ratio of partially filled grid box height to its surrounding neighbours", margin_coeff,  margin_coeff_set); CHKERRQ(ierr);
    if (!margin_coeff_set) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: Please specify scale coefficient for scaleMargin method.\n");
      CHKERRQ(ierr);
      PISMEnd();
    }
  }
    
  // calculate Href as average height of the neighbouring grounded grid cells
  // get mean ice thickness over adjacent grounded ice shelf neighbors

  // if a direct neighbour is grounded ice, then add H to H_average and increase N.
  if (m.grounded_ice(M.e)) { H_average += H.e; N++; }
  if (m.grounded_ice(M.w)) { H_average += H.w; N++; }
  if (m.grounded_ice(M.n)) { H_average += H.n; N++; }
  if (m.grounded_ice(M.s)) { H_average += H.s; N++; }

  if (N == 0) {
    // assume constant H_average for isolated grounded cells
    H_average = 200.0;
    vTestVar(i,j)=1;
//     SETERRQ(grid.com, 1, "N == 0;  call this only if a neighbor is grounded!\n");
  }
  // average by division by grounded neighbours
  H_average = H_average / N;
  if (HrefCalcMethod == "Hav"){
//     ierr = verbPrintf(2, grid.com,"Method Hav, H.s=%f,H.n=%f,H_av=%f at i,j=%d,%d\n",H.s,H.n,H_average,i,j); CHKERRQ(ierr);
    return H_average;
    
  }
  if (HrefCalcMethod == "scaleMargin"){
    // this method uses the average height of surrounding neighbours to calculate the height of the partially filled box and scales it down with coefficient margin_coeff

    // we here calculate Href from bed elevation, this is used for upward sloping beds
    ierr = vh.begin_access(); CHKERRQ(ierr);
    PetscReal h_average = 0.0;
    PetscReal H_average_FromBed = 0.0;
    N = 0;
    if (m.grounded_ice(M.e)) { h_average += vh(i+1,j); N++; }
    if (m.grounded_ice(M.w)) { h_average += vh(i-1,j); N++; }
    if (m.grounded_ice(M.n)) { h_average += vh(i,j+1); N++; }
    if (m.grounded_ice(M.s)) { h_average += vh(i,j-1); N++; }

    if (N == 0) {
    // assume constant H_average for isolated grounded cells
    h_average = 200.0;
    }

    // ice thickness is difference between elevation and bed
    H_average_FromBed = h_average / N - bed_ij;
    ierr = vh.end_access(); CHKERRQ(ierr);
    // decide on which Href to use, this is usually H_average for downward sloping and H_average_FromBed for upward sloping ground.

//     if( i==5 ){
//       ierr = verbPrintf(2, grid.com,"marginScale, H_average=%f, H_av_FromBed=%f at i,j=%d,%d\n",H_average, H_average_FromBed,i,j); CHKERRQ(ierr);
//     }
    H_average = PetscMin(H_average, H_average_FromBed);
//     if( i==5 ){
//       ierr = verbPrintf(2, grid.com,"marginScale, H_average chosen=%f at i,j=%d,%d\n",H_average, i,j); CHKERRQ(ierr);
//     }

    ierr = vTestVar.end_access(); CHKERRQ(ierr);
//     ierr = vh.end_access(); CHKERRQ(ierr);
    return H_average * margin_coeff;
        
  }


  if (HrefCalcMethod == "hav") {
    // this method is for testing, it is problematic for downward sloping beds
    // as it imposes a horizontal ice surface.
    ierr = vh.begin_access(); CHKERRQ(ierr);
    PetscReal h_average = 0.0;
    N = 0;
    if (m.grounded_ice(M.e)) { h_average += vh(i+1,j); N++; }
    if (m.grounded_ice(M.w)) { h_average += vh(i-1,j); N++; }
    if (m.grounded_ice(M.n)) { h_average += vh(i,j+1); N++; }
    if (m.grounded_ice(M.s)) { h_average += vh(i,j-1); N++; }

    if (N == 0) {
      SETERRQ(1, "N == 0;  call this only if a neighbor is grounded!\n");
    }
//     if( i==1 ){
//       ierr = verbPrintf(2, grid.com,"Method hav, h_average=%f,bed_ij=%f, H_average from Hav = %f, N=%d at i,j=%d,%d\n",h_average,bed_ij,H_average,N,i,j); CHKERRQ(ierr);
//     }

    // ice thickness is difference between elevation and bed
    
    H_average = h_average / N - bed_ij;
//     if( i==1 ){
//       ierr = verbPrintf(2, grid.com,"new H_average=%f at i,j=%d,%d\n",H_average,i,j); CHKERRQ(ierr);
//     }
    ierr = vh.end_access(); CHKERRQ(ierr);
    return H_average;
  }

  // shouldnt arrive here.
  SETERRQ(2, "no part_grid_ground method was chosen, check iMpartgridground.cc");
  

}

// //! \brief This is a copy of the Mask method in iMgeometry, merge the methods later.

// PetscErrorCode IceModel::update_mask_forpartgrid() {
//   PetscErrorCode ierr;
// 
//   if (ocean == PETSC_NULL) {  SETERRQ(grid.com, 1, "PISM ERROR: ocean == PETSC_NULL");  }
//   PetscReal sea_level;
//   ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
// 
//   GeometryCalculator gc(sea_level, *ice, config);
//   // construct the object mask with argument vMaskPartgrid, really needed here?
//   MaskQuery mask(vMaskPartGrid);
// 
//   ierr =    vH.begin_access();    CHKERRQ(ierr);
//   ierr =  vbed.begin_access();  CHKERRQ(ierr);
//   ierr = vMaskPartGrid.begin_access(); CHKERRQ(ierr);
// 
//   PetscInt GHOSTS = 2;
//   for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
//     for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
//       vMaskPartGrid(i, j) = gc.mask(vbed(i, j), vH(i,j));
//     } // inner for loop (j)
//   } // outer for loop (i)
// 
//   ierr =         vH.end_access(); CHKERRQ(ierr);
//   ierr =       vbed.end_access(); CHKERRQ(ierr);
//   ierr =      vMaskPartGrid.end_access(); CHKERRQ(ierr);
// 
//   return 0;
// }


//! \brief This method implements melting at grounded ocean margins.
/*!
  TODO: Give some detailed description here.
*/
