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
#include "PISMTime.hh"


//! \file iMfractures.cc implementing calculation of fracture density with PIK options -fractures.


PetscErrorCode IceModel::calculateFractureDensity() {
  const PetscScalar   dx = grid.dx, dy = grid.dy, Mx = grid.Mx, My = grid.My;
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### calculateFractureDensity() start \n");    CHKERRQ(ierr);

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


/*!
This routine provides the user to apply a icefree rift into a shelfsystem
*/
PetscErrorCode IceModel::applyRift() {
  PetscErrorCode ierr;
  const PetscScalar   dx = grid.dx, dy = grid.dy, Mx = grid.Mx, My = grid.My;
  ierr = verbPrintf(4,grid.com,"######### applyRift is called \n");    CHKERRQ(ierr);


  // -rift[riftnumber,time1,time2]
  // option -rift asks for three parameters, namely the rift number in the code as well as the initial and finish time of rift application
  
  const PetscBool  DEFAULT_DO_RIFT = PETSC_FALSE;
  const PetscInt DEFAULT_RIFT_NUMBER = 0;
  const PetscScalar DEFAULT_RIFT_TIME1 =0.0;
  const PetscScalar DEFAULT_RIFT_TIME2 =0.0;
  PetscScalar current_time = grid.time->current()/secpera;
  PetscInt rift_number = DEFAULT_RIFT_NUMBER;
  PetscScalar rift_time1 = DEFAULT_RIFT_TIME1;
  PetscScalar rift_time2 = DEFAULT_RIFT_TIME2;
  PetscInt    NparamRift=3;
  PetscReal   inarrayRift[3] = {DEFAULT_RIFT_NUMBER, DEFAULT_RIFT_TIME1,DEFAULT_RIFT_TIME2};
  PetscBool  RiftSet;

  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-rift", inarrayRift, &NparamRift, &RiftSet); CHKERRQ(ierr);
  if (RiftSet == PETSC_TRUE){
    rift_number = inarrayRift[0];
    rift_time1 = inarrayRift[1];
    rift_time2 = inarrayRift[2];
  }

  if (current_time>rift_time1 && current_time<rift_time2) {
    ierr = vH.begin_access(); CHKERRQ(ierr);

    //PetscScalar **Riftmask; 
    //ierr = vRiftmask.get_array(Riftmask); CHKERRQ(ierr);
    PetscBool RiftHere = PETSC_FALSE;

	
    for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) { 
      for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
	RiftHere=PETSC_FALSE;
	//Test 50km:
	if(rift_number==1){
	  if(i==floor(95*Mx/101) && j<floor(80*My/101)){ 
	    RiftHere= PETSC_TRUE;}}

       //Larsen 2km:
       //FIXME: define rift position more generally in terms of grid dimensions in order to use it for different resolutions
       if(rift_number==201) {//larsenA north and south
	  if((j<0.6*j+125 && i>0.6*i+122 && j>38)||( i<189 && i>186  && j>14)) { 
	    RiftHere= PETSC_TRUE; }}
       if(rift_number==202){//larsenB sorth north
	  if( i<0.8*j+85 && i>0.8*j+82 && j>60){  
	    RiftHere= PETSC_TRUE;}}
       if(rift_number==203){//larsenB north longer
	  if((i<2.0*j+14 && i>2.0*j+9 && j<62 && j>44)|| (i<0.8*j+85 && i>0.8*j+82 && j>60)) { 
	    RiftHere= PETSC_TRUE;}}
       if(rift_number==204){//larsenB south short
	  if(i<1.7*j-147 && i>1.7*j-151 && j>116){ 
	    RiftHere= PETSC_TRUE;}}	
       if(rift_number==205){//larsenB south longer
	  if((i<19*j-2130 && i>19*j-2165 && i>28 && i<51)||(i<1.7*j-147 && i>1.7*j-151 && j>116)) { 
	    RiftHere= PETSC_TRUE;}}
       if(rift_number==206){//larsenA south long
	  if(i<1.0*j+99 && i>1.0*j+94 && j>38){ 
	    RiftHere= PETSC_TRUE;}}
       if(rift_number==207){//larsenB north long
	  if(i<0.8*j+87 && i>0.8*j+83 && j>42){ 
	    RiftHere= PETSC_TRUE;}}
       if(rift_number==208){//larsenB south long
	  if(i<1.3*j-91 && i>1.3*j-95 && j>96){ 
	    RiftHere= PETSC_TRUE;}}
       if(rift_number==209){//larsenB north and south long
	  if((i<1.3*j-91 && i>1.3*j-95 && j>96)||(i<0.8*j+87 && i>0.8*j+83 && j>40)){ 
	    RiftHere= PETSC_TRUE;}}
      

       //ross 6,1km:
       if(rift_number==301) {//ross right up to byrd
	  if(i<2.684*j-283 && i>2.684*j-288 && i<85) { 
	    RiftHere= PETSC_TRUE; }}
       if(rift_number==302) {//ross right of ice rise
	  if(j==39 && i<90) { 
	    RiftHere= PETSC_TRUE; }}


       if(RiftHere)
          vH(i,j)=0.0;       
      }
    }
 

    RiftIsCut = PETSC_TRUE; // this ensures that a rift is cut ONLY on one single timestep
    ierr = verbPrintf(1,grid.com,"PISMPIK_INFO: Rift number %d is cut at year %f between %f and %f\n",rift_number,current_time,rift_time1,rift_time2);
    ierr = vH.end_access(); CHKERRQ(ierr); 
  }

  return 0;
}


