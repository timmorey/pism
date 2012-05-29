// Copyright (C) 2004--2012 Torsten Albrecht
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
#include <petscdmda.h>
#include "iceModel.hh"
#include "pism_signal.h"
#include "Mask.hh"
#include "PISMStressBalance.hh"


//! \file iMfractures.cc implementing calculation of fracture density with PIK options -fractures.


PetscErrorCode IceModel::calculateFractureDensity() {
  const PetscScalar   dx = grid.dx, dy = grid.dy, Mx = grid.Mx, My = grid.My;
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com, "######### calculateFractureDensity() start \n");    CHKERRQ(ierr);

  IceModelVec2S vFDnew = vWork2d[0];
  IceModelVec2S vFAnew = vWork2d[1];
    
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vFD.begin_access(); CHKERRQ(ierr);
  ierr = vFD.copy_to(vFDnew); CHKERRQ(ierr);
  ierr = vFDnew.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access();  CHKERRQ(ierr);
   
  ierr = stress_balance->get_principal_strain_rates(vPrinStrain1, vPrinStrain2); CHKERRQ(ierr);
  ierr = stress_balance->get_2D_stresses(txx, tyy, txy); CHKERRQ(ierr);
   
  ierr = txx.beginGhostComm(); CHKERRQ(ierr);
  ierr = txx.endGhostComm(); CHKERRQ(ierr);
  ierr = tyy.beginGhostComm(); CHKERRQ(ierr);
  ierr = tyy.endGhostComm(); CHKERRQ(ierr);
  ierr = txy.beginGhostComm(); CHKERRQ(ierr);
  ierr = txy.endGhostComm(); CHKERRQ(ierr);
    
  ierr = vPrinStrain1.beginGhostComm(); CHKERRQ(ierr);
  ierr = vPrinStrain1.endGhostComm(); CHKERRQ(ierr);
  ierr = vPrinStrain2.beginGhostComm(); CHKERRQ(ierr);
  ierr = vPrinStrain2.endGhostComm(); CHKERRQ(ierr);
   
  ierr = vPrinStrain1.begin_access(); CHKERRQ(ierr);
  ierr = vPrinStrain2.begin_access(); CHKERRQ(ierr);
  ierr = txx.begin_access(); CHKERRQ(ierr);
  ierr = tyy.begin_access(); CHKERRQ(ierr);
  ierr = txy.begin_access(); CHKERRQ(ierr);
  
  
  const bool dirichlet_bc = config.get_flag("ssa_dirichlet_bc");
  if (dirichlet_bc) {
    ierr = vBCMask.begin_access();  CHKERRQ(ierr);
    ierr = vBCvel.begin_access();  CHKERRQ(ierr);
  }
  
  const bool write_fd = config.get_flag("write_fd_fields");
  if (write_fd) {
    //double n_glen  = ice.exponent();
    ierr = vFG.begin_access();  CHKERRQ(ierr);
    ierr = vFH.begin_access();  CHKERRQ(ierr);
    ierr = vFE.begin_access();  CHKERRQ(ierr);
    ierr = vFA.begin_access();  CHKERRQ(ierr);
    ierr = vFA.copy_to(vFAnew); CHKERRQ(ierr);
    ierr = vFAnew.begin_access(); CHKERRQ(ierr);
   }
   
   IceModelVec2V *vel_advective;
   ierr = stress_balance->get_advective_2d_velocity(vel_advective); CHKERRQ(ierr);
   //IceModelVec2V vel = *vel_advective; // just an alias
   ierr = vel_advective->begin_access(); CHKERRQ(ierr);
  
   const bool do_part_grid = config.get_flag("part_grid");
   if  (do_part_grid == 0)
     ierr = verbPrintf(2, grid.com,"!!! PISM_WARNING: Fracture density only with -part_grid\n"); CHKERRQ(ierr);
  
   MaskQuery M(vMask);
   PetscScalar tempFD;


   /////////////////////////////////////////////////////////
   //get four options for calculation of fracture density.
   //1st: fracture growth constant gamma
   //2nd: fracture initiation stress threshold sigma_cr 
   //3rd: healing rate constant gamma_h
   //4th: healing starin rate threshold
   //more: T. Albrecht, A. Levermann; Fracture field for large-scale ice dynamics; (2012), 
   //Journal of Glaciology, Vol. 58, No. 207, 165-176, DOI: 10.3189/2012JoG11J191.

   PetscInt  Nparamf=4;
   PetscReal inarrayf[4] = {1.0, 7.0e4, 0.0, 2.0e-10};
  
   PetscBool  fractures_set;
   ierr = PetscOptionsGetRealArray(PETSC_NULL, "-fractures", inarrayf, &Nparamf, &fractures_set);
   CHKERRQ(ierr);
  
   if ((Nparamf > 4) || (Nparamf < 4)) {
     ierr = verbPrintf(1, grid.com,
                           "PISM ERROR: option -fractures provided with more or fewer than 4\n"
                           "            arguments ... ENDING ...\n");CHKERRQ(ierr);
     PISMEnd();
   }
   PetscReal gamma = inarrayf[0], initThreshold = inarrayf[1], gammaheal= inarrayf[2], healThreshold = inarrayf[3];
  
   ierr = verbPrintf(4, grid.com,"PISM-PIK INFO: fracture density is found with gamma=%.2f, sigma_cr=%.2f, gammah=%.2f and healing_cr=%e \n", gamma,initThreshold,gammaheal,healThreshold); CHKERRQ(ierr);
     
   
   PetscScalar fdBoundaryValue = 0.0;
   ierr = PetscOptionsGetScalar(PETSC_NULL, "-phi0", &fdBoundaryValue, PETSC_NULL);
     


   /////////////////////////////////////////////////////////////////////
   //get three options for flow enhancement due to fractures.
   //assume linear response function: E_fr = soft_rate * fracture_density + soft_offset if fracture_density > phi_init
   //set option as -fracture_softening 0.0,4.5,1.0
   PetscInt    Nparam=3;
   PetscReal   inarray[3] = {0.0,0.0,1.0};//phi_init, soft_rate, soft_offset
  
   PetscBool sd_set;
   ierr = PetscOptionsGetRealArray(PETSC_NULL, "-fracture_softening", inarray, &Nparam, &sd_set);
   CHKERRQ(ierr);
  
   PetscReal phi_init= inarray[0], soft_rate = inarray[1], soft_offset = inarray[2];
  
   if (sd_set) {
     if ((Nparam > 3) || (Nparam < 3)) {
        ierr = verbPrintf(1, grid.com,
                           "PISM ERROR: option -fracture_softening provided with more or fewer than 3\n"
                           "            arguments ... ENDING ...\n");CHKERRQ(ierr);
        PISMEnd();
     }
     ierr = verbPrintf(4, grid.com,"PISM-PIK INFO: fracture_softening mode is set with phi_min=%.2f, ms=%.2f and ns=%.2f\n", phi_init,soft_rate,soft_offset); CHKERRQ(ierr);
  
   }
   
    
   for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
     for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

       tempFD=0;
       //SSA: v . grad memField
       tempFD +=  (*vel_advective)(i,j).u * ((*vel_advective)(i,j).u<0 ? vFD(i+1,j)-vFD(i, j):vFD(i, j)-vFD(i-1, j))/dx; 
       tempFD +=  (*vel_advective)(i,j).v * ((*vel_advective)(i,j).v<0 ? vFD(i,j+1)-vFD(i, j):vFD(i, j)-vFD(i, j-1))/dy; 

       vFDnew(i,j)-= tempFD * dt;
       
       //source
       //van mises
       PetscScalar T1 =0.5*(txx(i,j)+tyy(i,j))+sqrt(0.25*PetscSqr(txx(i,j)-tyy(i,j))+PetscSqr(txy(i,j)));//Pa
       PetscScalar T2 =0.5*(txx(i,j)+tyy(i,j))-sqrt(0.25*PetscSqr(txx(i,j)-tyy(i,j))+PetscSqr(txy(i,j)));//Pa
       PetscScalar sigmat=sqrt(PetscSqr(T1)+PetscSqr(T2)-T1*T2);

  
       if (sigmat > initThreshold) 
          vFDnew(i,j)+= gamma*(vPrinStrain1(i,j)-0.0)*dt*(1-vFDnew(i,j));       
       //if (sigmat > initThreshold && vPrinStrain1(i,j)>3.0e-10) 
          //vFDnew(i,j)+= gamma*(vPrinStrain1(i,j)-3.0e-10)*dt*(1-vFDnew(i,j));
       //if (vPrinStrain1(i,j)>initThreshold) 
          //vFDnew(i,j)+= gamma*(vPrinStrain1(i,j))*dt*(1-vFDnew(i,j)));
           
       //healing
       if (vPrinStrain1(i,j) < healThreshold) {
           vFDnew(i,j)+= gammaheal*(vPrinStrain1(i,j)-healThreshold)*dt;//*(1-vFD(i,j)); 
       }
       //vFDnew(i,j)+= gammaheal*(-healThreshold)*dt;
       
       // write related fracture quantities to nc-file
       // if option -write_fd_fields is set
       if (write_fd && vH(i,j)>0.0) {
         // fracture growth rate
         if (sigmat > initThreshold) {
           vFG(i,j)=gamma*(vPrinStrain1(i,j)-0.0)*(1-vFDnew(i,j)); }
         else {
           vFG(i,j)=0.0;}
         
         // fracture healing rate      
         if (vPrinStrain1(i,j) < healThreshold){
           vFH(i,j)=gammaheal*(vPrinStrain1(i,j)-healThreshold); }
         else {
           vFH(i,j)=0.0;}
         //vFH(i,j)=gammaheal*(-healThreshold);//}              
    
         //fracture age since fracturing occured
         if (sigmat <= initThreshold){// && (vFA(i+1,j)>0.0 || vFA(i-1,j)>0.0 || vFA(i,j+1)>0.0 || vFA(i,j-1)>0.0)) {
           vFAnew(i,j) -= dt * (*vel_advective)(i,j).u * ((*vel_advective)(i,j).u<0 ? vFA(i+1,j)-vFA(i, j):vFA(i, j)-vFA(i-1, j))/dx; 
           vFAnew(i,j) -= dt * (*vel_advective)(i,j).v * ((*vel_advective)(i,j).v<0 ? vFA(i,j+1)-vFA(i, j):vFA(i, j)-vFA(i, j-1))/dy;  
         }
         if (sigmat <= initThreshold)// && vFAnew(i,j)>0.0)
           vFAnew(i,j)+= dt/secpera;
         else if (sigmat > initThreshold) 
           vFAnew(i,j) = 0.0;

         // additional flow enhancement due to fracture softening
         PetscScalar softening = soft_offset + soft_rate*vFDnew(i,j);        
         if (vH(i,j)>0.0)
           vFE(i,j)=1.0/pow(softening,1/3.0);
         else
           vFE(i,j)=0.0;
       }  
   
       //boundary condition
       if (dirichlet_bc) {
         if (vBCMask.as_int(i,j) == 1 && (vBCvel(i,j).u != 0.0 || vBCvel(i,j).v != 0.0))
            vFDnew(i,j)=fdBoundaryValue;
       }
             
      //bounding        
      if (vFDnew(i,j)<0.0)
        vFDnew(i,j)=0.0;
      if (vFDnew(i,j)>1.0)
         vFDnew(i,j)=1.0;
             
      if (vH(i,j)==0.0)
        vFDnew(i,j)=0.0;
 
     }
   } 
    
    ierr = vel_advective->end_access(); CHKERRQ(ierr);
    ierr = vH.end_access(); CHKERRQ(ierr);
    ierr = vFD.end_access(); CHKERRQ(ierr);
    ierr = vFDnew.end_access(); CHKERRQ(ierr);
    ierr = vMask.end_access();  CHKERRQ(ierr);
  
   if (dirichlet_bc) {
     ierr = vBCMask.end_access();  CHKERRQ(ierr);
     ierr = vBCvel.end_access();  CHKERRQ(ierr);
    }
    
   if (write_fd) {
     ierr = vFG.end_access();  CHKERRQ(ierr);
     ierr = vFH.end_access();  CHKERRQ(ierr);
     ierr = vFE.end_access();  CHKERRQ(ierr);
     ierr = vFA.end_access();  CHKERRQ(ierr);
     ierr = vFAnew.end_access(); CHKERRQ(ierr);
     ierr = vFAnew.beginGhostComm(vFA); CHKERRQ(ierr);
     ierr = vFAnew.endGhostComm(vFA); CHKERRQ(ierr);
   }
     
   ierr = vPrinStrain1.end_access(); CHKERRQ(ierr);
   ierr = vPrinStrain2.end_access(); CHKERRQ(ierr);
   ierr = txx.end_access(); CHKERRQ(ierr);
   ierr = tyy.end_access(); CHKERRQ(ierr);
   ierr = txy.end_access(); CHKERRQ(ierr);
   
   ierr = vFDnew.beginGhostComm(vFD); CHKERRQ(ierr);
   ierr = vFDnew.endGhostComm(vFD); CHKERRQ(ierr);
 
  return 0;
}
