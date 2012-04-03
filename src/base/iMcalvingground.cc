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

//! \brief This method implements melting at grounded ocean margins.
/*!
  TODO: Give some detailed description here.
*/

PetscErrorCode IceModel::groundedCalving() {
  PetscErrorCode ierr;
  bool calvMethod_set;
  string GroundedCalvingMethod;
  set<string> GroundCalv_choices;
  GroundCalv_choices.insert("constant");
  GroundCalv_choices.insert("eigen");
  GroundCalv_choices.insert("eigen_and_const");
  ierr = PISMOptionsList(grid.com,"-grounded_calving", "specifies the grounded calving calculation method",
         GroundCalv_choices, "constant", GroundedCalvingMethod, calvMethod_set); CHKERRQ(ierr);
  if (!calvMethod_set) {
    ierr = PetscPrintf(grid.com,
           "PISM ERROR: Please specify an method for grounded calving.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }
  
  if (GroundedCalvingMethod == "constant"){
    ierr = groundedCalvingConst(); CHKERRQ(ierr);
  }
  if (GroundedCalvingMethod == "eigen"){
    ierr = groundedEigenCalving(); CHKERRQ(ierr);
  }
  if (GroundedCalvingMethod == "eigen_and_const"){
    ierr = groundedEigenCalving(); CHKERRQ(ierr);
    ierr = groundedCalvingConst(); CHKERRQ(ierr);    
  }
  
  return 0;
}


//! \brief This method implements melting at grounded ocean margins.
/*!
  TODO: Give some detailed description here.
*/
PetscErrorCode IceModel::groundedEigenCalving() {
  const PetscScalar   dx = grid.dx, dy = grid.dy;
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### groundedEigenCalving() start \n");    CHKERRQ(ierr);

  PetscScalar my_eigencalv_ground_flux = 0.0;
  bool eigcalv_ground_factor_set;
  PetscReal eigcalv_ground_factor;
  ierr = PISMOptionsReal("-eigcalv_ground_factor", "specifies eigen calving factor for grounded margins.", eigcalv_ground_factor,  eigcalv_ground_factor_set); CHKERRQ(ierr);
  if (!eigcalv_ground_factor_set) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Please specify eigen calving factor for grounded margins.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }
  bool thresh_coeff_set;
  PetscReal thresh_coeff = 1.0;
  ierr = PISMOptionsReal("-thresh_coeff", "specifies a coefficient to avoid oscillations between HrefG and full cell", thresh_coeff, thresh_coeff_set); CHKERRQ(ierr);

  bool landeigencalving;
  ierr = PISMOptionsIsSet("-landeigencalving", "Use eigenCalvingGround also on land above SL.", landeigencalving); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density");
  double ice_rho   = config.get("ice_density");
  double rhofrac   = ice_rho/ocean_rho;
  // Distance (grid cells) from calving front where strain rate is evaluated
  PetscInt offset = 2;

  PetscReal sea_level = 0;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(2, "PISM ERROR: ocean == NULL"); }

  if (PetscAbs(dx - dy)/PetscMin(dx,dy) > 1e-2) {
    ierr = PetscPrintf(grid.com,
      "PISMPIK_ERROR: -eigen_calving using a non-square grid cell does not work (yet);\n"
                       "  since it has no direction!!!\n, dx = %f, dy = %f, rel. diff = %f",dx,dy,PetscAbs(dx - dy)/PetscMax(dx,dy));
    PISMEnd();
  }

  MaskQuery mask(vMask);

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vHavgGround.begin_access(); CHKERRQ(ierr);
  ierr = vHrefGround.begin_access(); CHKERRQ(ierr);
  ierr = vJustGotFullCell.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = vPrinStrain1.begin_access(); CHKERRQ(ierr);
  ierr = vPrinStrain2.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  IceModelVec2S vDiffCalvHeight = vWork2d[1];
  ierr = vDiffCalvHeight.set(0.0); CHKERRQ(ierr);
  ierr = vDiffCalvHeight.begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {   
      // we should substitute these definitions, can we have a more flexible mask class
      // that we can use to create a mask object that is not the default mask and can
      // handle part grid cells?
      bool grounded_ice     = vH(i,j) > 0.0 && vbed(i,j) > (sea_level - rhofrac*vH(i,j));
      bool part_grid_cell   = vHrefGround(i,j) > 0.0;
      bool below_sealevel   = (vbed(i,j) + sea_level) < 0;
      bool free_ocean       = vH(i,j) == 0.0 && below_sealevel && !part_grid_cell;
      bool grounded_e = vH(i+1,j)>0.0 && vbed(i+1,j)>(sea_level-rhofrac*vH(i+1,j)) && (vbed(i+1,j)+sea_level)<0;
      bool grounded_w = vH(i-1,j)>0.0 && vbed(i-1,j)>(sea_level-rhofrac*vH(i-1,j)) && (vbed(i-1,j)+sea_level)<0;
      bool grounded_n = vH(i,j+1)>0.0 && vbed(i,j+1)>(sea_level-rhofrac*vH(i,j+1)) && (vbed(i,j+1)+sea_level)<0;
      bool grounded_s = vH(i,j-1)>0.0 && vbed(i,j-1)>(sea_level-rhofrac*vH(i,j-1)) && (vbed(i,j-1)+sea_level)<0;
      bool next_to_grounded = grounded_e || grounded_w || grounded_n || grounded_s;
      // ocean front is where no partially filled grid cell is in front.
      bool at_ocean_front_e = ( vH(i+1,j) == 0.0 && ((vbed(i+1,j) + sea_level) < 0 && vHrefGround(i+1,j) == 0.0) );
      bool at_ocean_front_w = ( vH(i-1,j) == 0.0 && ((vbed(i-1,j) + sea_level) < 0 && vHrefGround(i-1,j) == 0.0) );
      bool at_ocean_front_n = ( vH(i,j+1) == 0.0 && ((vbed(i,j+1) + sea_level) < 0 && vHrefGround(i,j+1) == 0.0) );
      bool at_ocean_front_s = ( vH(i,j-1) == 0.0 && ((vbed(i,j-1) + sea_level) < 0 && vHrefGround(i,j-1) == 0.0) );
      bool at_ocean_front   = at_ocean_front_e || at_ocean_front_w || at_ocean_front_n || at_ocean_front_s;

                               
      if( part_grid_cell && (at_ocean_front && below_sealevel || landeigencalving)){
        PetscScalar dHref = 0.0, Face = 0.0, calvrateHorizontal = 0.0,
        eigen1 = 0.0, eigen2 = 0.0, 
        eigenCalvOffset = 0.0; // if it's not exactly the zero line of
                               // transition from compressive to extensive flow
                               // regime;
        // Counting adjacent grounded boxes (with distance "offset")
        PetscInt M = 0, N = 0;
        
        Face = 1.0/dx;
        // is this a good idea for non quadratic grid handling?
//         if ( at_ocean_front_e ) { Face+= 1.0/dx; N++;}
//         if ( at_ocean_front_w ) { Face+= 1.0/dx; N++;}
//         if ( at_ocean_front_n ) { Face+= 1.0/dy; N++;}
//         if ( at_ocean_front_s ) { Face+= 1.0/dy; N++;}
//         if (N > 0) Face = Face/N;
        
        // make this less rough if in use for future.
        if ( landeigencalving ) Face = 1.0/dx;
                
        if ( mask.grounded_ice(i + offset, j) && !mask.ice_margin(i + offset, j)){
          eigen1 += vPrinStrain1(i + offset, j);
          eigen2 += vPrinStrain2(i + offset, j);
          M += 1;
        }
        if ( mask.grounded_ice(i - offset, j) && !mask.ice_margin(i - offset, j)){
          eigen1 += vPrinStrain1(i - offset, j);
          eigen2 += vPrinStrain2(i - offset, j);
          M += 1;
        }
        if ( mask.grounded_ice(i, j + offset) && !mask.ice_margin(i , j + offset)){
          eigen1 += vPrinStrain1(i, j + offset);
          eigen2 += vPrinStrain2(i, j + offset);
          M += 1;
        }
        if ( mask.grounded_ice(i, j - offset) && !mask.ice_margin(i , j - offset)){
          eigen1 += vPrinStrain1(i, j - offset);
          eigen2 += vPrinStrain2(i, j - offset);
          M += 1;
        }
        if (M > 0) {
          eigen1 /= M;
          eigen2 /= M;
        }


        // calving law
        if ( eigen2 > eigenCalvOffset && eigen1 > 0.0) { // if spreading in all directions
          calvrateHorizontal = eigcalv_ground_factor * eigen1 * (eigen2 - eigenCalvOffset);
          // eigen1 * eigen2 has units [s^ - 2] and calvrateHorizontal [m*s^1]
          // hence, eigcalv_ground_factor has units [m*s]
        } else calvrateHorizontal = 0.0;

        // calculate mass loss with respect to the associated ice thickness and the grid size:
        PetscScalar calvrate = calvrateHorizontal * vHavgGround(i,j) * Face; // in m/s
        // dHref corresponds to the height we have to cut off to mimic a
        // a constant horizontal retreat of a part grid cell.
        // volume_partgrid = Href * dx*dy
        // area_partgrid   = volume_partgrid/Havg = Href/Havg * dx*dy
        // calv_velocity   = const * d/dt(area_partgrid/dy) = const * dHref/dt * dx/Havg    
        dHref = calvrate * dt;

        // all further calving originates from dHref, so this line is enough.
        my_eigencalv_ground_flux -= dHref;
        if( vHrefGround(i,j) > dHref ){
          // enough ice to calv from partial cell
          vHrefGround(i,j) -= dHref;
        } else {
        PetscInt N = 0;
        // count grounded neighbours
        if ( grounded_e && vJustGotFullCell(i+1,j) != 1. ) N++;
        if ( grounded_w && vJustGotFullCell(i-1,j) != 1. ) N++;
        if ( grounded_n && vJustGotFullCell(i,j+1) != 1. ) N++;
        if ( grounded_s && vJustGotFullCell(i,j-1) != 1. ) N++;
        // kill partial cell and save for redistribution to grounded neighbours
        // isolated (N==0) PGG cell at ocean front is killed without redistribution.
        if (N > 0) vDiffCalvHeight(i,j) = (dHref - vHrefGround(i,j))/N;
        vHrefGround(i,j)     = 0.0;
        }
      }
    }
  }

  ierr = vDiffCalvHeight.end_access(); CHKERRQ(ierr);

  ierr = vDiffCalvHeight.beginGhostComm(); CHKERRQ(ierr);
  ierr = vDiffCalvHeight.endGhostComm(); CHKERRQ(ierr);

  ierr = vDiffCalvHeight.begin_access(); CHKERRQ(ierr);
  ierr = vHrefThresh.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      PetscScalar restCalvHeight = 0.0;
      bool grounded_ice   = vH(i,j) > 0.0 && vbed(i,j) > (sea_level - rhofrac*vH(i,j));
      bool below_sealevel = (vbed(i,j) + sea_level) < 0;
      // if grounded and DiffCalv in a neighbouring cell
      if ( grounded_ice && below_sealevel && vJustGotFullCell(i,j) != 1. &&
           (vDiffCalvHeight(i + 1, j) > 0.0 || vDiffCalvHeight(i - 1, j) > 0.0 ||
            vDiffCalvHeight(i, j + 1) > 0.0 || vDiffCalvHeight(i, j - 1) > 0.0 )) {

        restCalvHeight = vDiffCalvHeight(i + 1, j) + vDiffCalvHeight(i - 1, j) +
                         vDiffCalvHeight(i, j + 1) + vDiffCalvHeight(i, j - 1);

        vHrefGround(i, j) = vH(i, j) - restCalvHeight; // in m
        PetscSynchronizedPrintf(grid.com,"make Hnew=%e a Href=%e cell with rCalv= %e at i=%d, j=%d\n",vH(i, j),vHrefGround(i, j), restCalvHeight,i,j);

        vHrefThresh(i,j) = vH(i, j) * thresh_coeff;
        vHnew(i, j)      = 0.0;

        if(vHrefGround(i, j) < 0.0) { // i.e. terminal floating ice grid cell has calved off completely.
          // We do not account for further calving ice-inwards!
          vHrefGround(i, j) = 0.0;
        }
      }
    }
  }

  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  ierr = vHavgGround.end_access(); CHKERRQ(ierr);
  ierr = vHrefGround.end_access(); CHKERRQ(ierr);
  ierr = vHrefThresh.end_access(); CHKERRQ(ierr);
  ierr = vDiffCalvHeight.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vJustGotFullCell.end_access(); CHKERRQ(ierr);
  ierr = vPrinStrain1.end_access(); CHKERRQ(ierr);
  ierr = vPrinStrain2.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&my_eigencalv_ground_flux, &eigencalv_ground_flux, grid.com); CHKERRQ(ierr);

  // FIXME: use corrected cell areas (when available)
  PetscScalar ice_density = config.get("ice_density"),
  factor = ice_density * (dx * dy) / dt;
  eigencalv_ground_flux *= factor;

  return 0;
}


//! \brief This method implements melting at grounded ocean margins.
/*!
  TODO: Give some detailed description here.
*/
PetscErrorCode IceModel::groundedCalvingConst() {
  const PetscScalar   dx = grid.dx, dy = grid.dy;
  PetscScalar my_thkcalv_ground_flux = 0.0;
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### groundedCalvingConst() start \n");    CHKERRQ(ierr);

  bool melt_factor_set;
  PetscReal ocean_melt_factor;
  ierr = PISMOptionsReal("-ocean_melt_factor", "specifies constant melt factor for oceanic melt at grounded margins.", ocean_melt_factor,  melt_factor_set); CHKERRQ(ierr);
  if (!melt_factor_set) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Please specify melt coefficient for ocean melt.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }
  bool thresh_coeff_set;
  PetscReal thresh_coeff = 1.0;
  ierr = PISMOptionsReal("-thresh_coeff", "specifies a coefficient to avoid oscillations between HrefG and full cell", thresh_coeff, thresh_coeff_set); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density");
  double ice_rho   = config.get("ice_density");
  double rhofrac   = ice_rho/ocean_rho;
  
  PetscReal sea_level = 0;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(2, "PISM ERROR: ocean == NULL"); }

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vHavgGround.begin_access(); CHKERRQ(ierr);
  ierr = vHrefGround.begin_access(); CHKERRQ(ierr);
  ierr = vJustGotFullCell.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  
  IceModelVec2S vDiffCalvHeight = vWork2d[1];
  ierr = vDiffCalvHeight.set(0.0); CHKERRQ(ierr);
  ierr = vDiffCalvHeight.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      // we should substitute these definitions, can we have a more flexible mask class
      // that we can use to create a mask object that is not the default mask and can
      // handle part grid cells?
      bool grounded_ice     = vH(i,j) > 0.0 && vbed(i,j) > (sea_level - rhofrac*vH(i,j));
      bool part_grid_cell   = vHrefGround(i,j) > 0.0;
      bool below_sealevel   = (vbed(i,j) + sea_level) < 0;
      bool free_ocean       = vH(i,j) == 0.0 && below_sealevel && !part_grid_cell;

      
      bool grounded_e = vH(i+1,j)>0.0 && vbed(i+1,j)>(sea_level-rhofrac*vH(i+1,j)) && (vbed(i+1,j)+sea_level)<0;
      bool grounded_w = vH(i-1,j)>0.0 && vbed(i-1,j)>(sea_level-rhofrac*vH(i-1,j)) && (vbed(i-1,j)+sea_level)<0;
      bool grounded_n = vH(i,j+1)>0.0 && vbed(i,j+1)>(sea_level-rhofrac*vH(i,j+1)) && (vbed(i,j+1)+sea_level)<0;
      bool grounded_s = vH(i,j-1)>0.0 && vbed(i,j-1)>(sea_level-rhofrac*vH(i,j-1)) && (vbed(i,j-1)+sea_level)<0;
      bool next_to_grounded = grounded_e || grounded_w || grounded_n || grounded_s;
      // ocean front is where no partially filled grid cell is in front.
      bool at_ocean_front_e = ( vH(i+1,j) == 0.0 && ((vbed(i+1,j) + sea_level) < 0 && vHrefGround(i+1,j) == 0.0) );
      bool at_ocean_front_w = ( vH(i-1,j) == 0.0 && ((vbed(i-1,j) + sea_level) < 0 && vHrefGround(i-1,j) == 0.0) );
      bool at_ocean_front_n = ( vH(i,j+1) == 0.0 && ((vbed(i,j+1) + sea_level) < 0 && vHrefGround(i,j+1) == 0.0) );
      bool at_ocean_front_s = ( vH(i,j-1) == 0.0 && ((vbed(i,j-1) + sea_level) < 0 && vHrefGround(i,j-1) == 0.0) );
      bool at_ocean_front   = at_ocean_front_e || at_ocean_front_w || at_ocean_front_n || at_ocean_front_s;

      if( part_grid_cell && at_ocean_front && below_sealevel ){
        PetscReal dHref = 0.0;
        if ( at_ocean_front_e ) { dHref+= 1/dx; }
        if ( at_ocean_front_w ) { dHref+= 1/dx; }
        if ( at_ocean_front_n ) { dHref+= 1/dy; }
        if ( at_ocean_front_s ) { dHref+= 1/dy; }
        // dHref corresponds to the height we have to cut off to mimic a
        // a constant horizontal retreat of a part grid cell.
        // volume_partgrid = Href * dx*dy
        // area_partgrid   = volume_partgrid/Havg = Href/Havg * dx*dy
        // calv_velocity   = const * d/dt(area_partgrid/dy) = const * dHref/dt * dx/Havg
        dHref = dHref * vHavgGround(i,j) * ocean_melt_factor * dt/secpera;
        // all further calving originates from dHref, so this line is enough.
        my_thkcalv_ground_flux -= dHref;        
        ierr = verbPrintf(2, grid.com,"dHref=%e at i=%d, j=%d\n",dHref,i,j); CHKERRQ(ierr);
        if( vHrefGround(i,j) > dHref ){
          // enough ice to calv from partial cell
          vHrefGround(i,j) -= dHref;
        } else {
        PetscInt N = 0;
        // count grounded neighbours
        if ( grounded_e && vJustGotFullCell(i+1,j) != 1. ) N++;
        if ( grounded_w && vJustGotFullCell(i-1,j) != 1. ) N++;
        if ( grounded_n && vJustGotFullCell(i,j+1) != 1. ) N++;
        if ( grounded_s && vJustGotFullCell(i,j-1) != 1. ) N++;
        // kill partial cell and save for redistribution to grounded neighbours
        // isolated (N==0) PGG cell at ocean front is killed without redistribution.
        if (N > 0) vDiffCalvHeight(i,j) = (dHref - vHrefGround(i,j))/N;
        vHrefGround(i,j)     = 0.0;
        }
      }
    }
  }

  ierr = vDiffCalvHeight.end_access(); CHKERRQ(ierr);

  ierr = vDiffCalvHeight.beginGhostComm(); CHKERRQ(ierr);
  ierr = vDiffCalvHeight.endGhostComm(); CHKERRQ(ierr);

  ierr = vDiffCalvHeight.begin_access(); CHKERRQ(ierr);
  ierr = vHrefThresh.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      PetscScalar restCalvHeight = 0.0;
      bool grounded_ice = vH(i,j) > 0.0 && vbed(i,j) > (sea_level - rhofrac*vH(i,j));
      bool below_sealevel   = (vbed(i,j) + sea_level) < 0;
      // if grounded and DiffCalv in a neighbouring cell
      if ( grounded_ice && below_sealevel && vJustGotFullCell(i,j) != 1. &&
           (vDiffCalvHeight(i + 1, j) > 0.0 || vDiffCalvHeight(i - 1, j) > 0.0 ||
            vDiffCalvHeight(i, j + 1) > 0.0 || vDiffCalvHeight(i, j - 1) > 0.0 )) {

        restCalvHeight = vDiffCalvHeight(i + 1, j) + vDiffCalvHeight(i - 1, j) +
                         vDiffCalvHeight(i, j + 1) + vDiffCalvHeight(i, j - 1);

        vHrefGround(i, j) = vH(i, j) - restCalvHeight; // in m
        PetscSynchronizedPrintf(grid.com,"make Hnew=%e a Href=%e cell with rCalv= %e at i=%d, j=%d\n",vH(i, j),vHrefGround(i, j), restCalvHeight,i,j);
      
        vHrefThresh(i,j) = vH(i, j) * thresh_coeff;
        vHnew(i, j)      = 0.0;

        if(vHrefGround(i, j) < 0.0) { // i.e. terminal floating ice grid cell has calved off completely.
          // We do not account for further calving ice-inwards!
          vHrefGround(i, j) = 0.0;
        }
      }
    }
  }

  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  ierr = vHavgGround.end_access(); CHKERRQ(ierr);
  ierr = vHrefGround.end_access(); CHKERRQ(ierr);
  ierr = vHrefThresh.end_access(); CHKERRQ(ierr);
  ierr = vDiffCalvHeight.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vJustGotFullCell.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&my_thkcalv_ground_flux, &thkcalv_ground_flux, grid.com); CHKERRQ(ierr);
  // FIXME: use corrected cell areas (when available)
  PetscScalar ice_density = config.get("ice_density"),
  factor = ice_density * (dx * dy) / dt;
  thkcalv_ground_flux *= factor;
  
  return 0;
}

