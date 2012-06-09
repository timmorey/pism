// Copyright (C) 2011, 2012 Constantine Khroulev
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

#include "POMeltingParam.hh"

PetscErrorCode POMeltingParam3eqn::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 3eqn melting parametrizatin ocean model\n"
                    "reading ocean temperature and salinity from a file.\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing",
                        "absolute temperature at ice shelf base",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                       "salinity at ice shelf base",
                       "psu", "salinity"); CHKERRQ(ierr);

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(1, "ERROR: ice thickness is not available"); }
  bed = dynamic_cast<IceModelVec2S*>(vars.get("topg"));
  if (!ice_thickness) { SETERRQ(1, "ERROR: topography is not available"); }

  return 0;
}

PetscErrorCode POMeltingParam3eqn::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr = update_internal(t_years, dt_years); CHKERRQ(ierr);

  if (enable_time_averaging) {
    ierr = mass_flux.average(t, dt); CHKERRQ(ierr);
    ierr = temp.average(t, dt); CHKERRQ(ierr);
  } else {
    ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
    ierr = temp.get_record_years(t); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode POMeltingParam3eqn::shelf_base_temperature(IceModelVec2S &result) {

  PetscErrorCode ierr;

  const PetscScalar T0 = config.get("water_melting_point_temperature"), // K
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME issue #15
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POMeltingParam3eqn::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

//   PetscReal L = config.get("water_latent_heat_fusion"),
    PetscReal rho_ocean = config.get("sea_water_density"),
              rho_ice   = config.get("ice_density"),
              meltrate_3eqn, meltrate_dummy, temp_insitu;

  PetscReal reference_pressure = 1.01325; // pressure of atmosphere in bar

  PetscScalar **H, **topgg;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = bed->get_array(topgg);CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux.begin_access(); CHKERRQ(ierr);
  


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // compute T_f[i][j] according to beckmann_goosse03, which has the
      // meaning of the freezing temperature of the ocean water directly
      // under the shelf, (of salinity 35psu) [this is related to the
      // Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for
      // details]
      PetscScalar shelfbaseelev = - (rho_ice / rho_ocean) * H[i][j];
      PetscReal thetao = temp(i,j) - 273.15; // in degC
      PetscReal sal  = mass_flux(i,j);
      // press in bar
      PetscReal press = rho_ocean * -1.*shelfbaseelev/1000. + reference_pressure;
//       ierr = verbPrintf(2, grid.com, "sal=%e, press=%e,temp=%e,\n",
//                     sal,press,temp(i,j)); CHKERRQ(ierr);

      potit(sal,thetao,press,reference_pressure, temp_insitu);
      PetscReal rhow = rho_ocean;
      PetscReal rhoi = rho_ice;
      PetscReal rho  = rho_ocean;
      PetscReal zice = -1*PetscAbs(shelfbaseelev);

//       ierr = verbPrintf(2, grid.com, "temp=%e, salt=%e\n",
//                         temp_insitu, sal); CHKERRQ(ierr);
      
      cavity_heat_water_fluxes_3eq(sal, temp_insitu, zice, rhow, rhoi, rho, meltrate_3eqn);

      result(i,j) = -1*meltrate_3eqn;

    }
  }

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POMeltingParam3eqn::cavity_heat_water_fluxes_3eq( PetscReal sal,
               PetscReal temp_insitu, PetscReal zice, PetscReal rhow,
               PetscReal rhoi, PetscReal rho, PetscReal &meltrate){

  // The three-equation model of ice-shelf ocean interaction (Hellmer and Olbers, 1989).
  // Code derived from BRIOS subroutine iceshelf (which goes back to H.Hellmer's 2D ice shelf model code)
  // and adjusted for use in FESOM by Ralph Timmermann, 16.02.2011
  // adapted for PISM by matthias.mengel@pik-potsdam.de
  
  PetscErrorCode ierr;
  PetscReal rhor, heat_flux, water_flux;
  PetscReal gats1, gats2, gas, gat;
  PetscReal ep1,ep2,ep3,ep4,ep5;
  PetscReal ex1,ex2,ex3,ex4,ex5;
  PetscReal vt1,sr1,sr2,sf1,sf2,tf1,tf2,tf,sf,seta,re;
  PetscInt n, n3, nk;

  PetscReal rp =   0.;                   //reference pressure
  PetscReal a   = -0.0575;                  //Foldvik&Kvinge (1974)
  PetscReal b   =  0.0901;
  PetscReal c   =  7.61e-4;

  PetscReal pr  =  13.8;                     //Prandtl number      [dimensionless]
  PetscReal sc  =  2432.;                  //Schmidt number      [dimensionless]
  PetscReal ak  =  2.50e-3;                 //dimensionless drag coeff.
  PetscReal sak1=  sqrt(ak);
  PetscReal un  =  1.95e-6;                //kinematic viscosity [m2/s]
  PetscReal pr1 =  pow(pr,(2./3.));              //Jenkins (1991)
  PetscReal sc1 =  pow(sc,(2./3.));  
  PetscReal tob=  -20.;                        //temperatur at the ice surface
//   PetscReal rhoi=  920.;                      //mean ice density
  PetscReal cpw =  4180.0;                    //Barnier et al. (1995)
  PetscReal lhf =  3.33e+5;                   //latent heat of fusion
  PetscReal atk =  273.15;                    //0 deg C in Kelvin
  //FIXME: can use PISMs surface temp for tob?
  PetscReal cpi =  152.5+7.122*(atk+tob);     //Paterson:"The Physics of Glaciers"

  PetscReal L    = 334000.;                   // [J/Kg]

// Prescribe the turbulent heat and salt transfer coeff. GAT and GAS

  gat  = 1.00e-4;   //[m/s] RG3417
  gas  = 5.05e-7;   //[m/s] RG3417

 // Calculate
 // density in the boundary layer: rhow
 // and interface pressure pg [dbar],
 // Solve a quadratic equation for the interface salinity sb
 // to determine the melting/freezing rate seta.

// FIXME: need to calculate water density instead of const value.
//  call fcn_density(thetao,sal,zice,rho)
//  rhow = density_0+rho  //was rhow= rho0+rho(i,j,N)
//  rhow = 1028.0;
 rhor= rhoi/rhow;

 ep1 = cpw*gat;
 ep2 = cpi*gas;
 ep3 = lhf*gas;
 ep4 = b+c*zice;
 ep5 = gas/rhor;
// negative heat flux term in the ice (due to -kappa/D)
 ex1 = a*(ep1-ep2);
 ex2 = ep1*(ep4-temp_insitu)+ep2*(tob+a*sal-ep4)-ep3;
 ex3 = sal*(ep2*(ep4-tob)+ep3);
 ex4 = ex2/ex1;
 ex5 = ex3/ex1;

 sr1 = 0.25*ex4*ex4-ex5;
 sr2 = -0.5*ex4;
 sf1 = sr2+sqrt(sr1);
 tf1 = a*sf1+ep4;
 sf2 = sr2-sqrt(sr1);
 tf2 = a*sf2+ep4;

// Salinities < 0 psu are not defined, therefore pick the positive of the two solutions:
 if(sf1 > 0.){
  tf = tf1;
  sf = sf1;
 }else{
  tf = tf2;
  sf = sf2;
 }

  // Calculate the melting/freezing rate [m/s]
  seta = ep5*(1.0-sal/sf);     //rt thinks this is not needed

  heat_flux  = rhow*cpw*gat*(temp_insitu-tf);      // [W/m2]
  water_flux =          gas*(sf-sal)/sf;   // [m/s]

  meltrate = seta;
//   ierr = verbPrintf(2, grid.com, "meltrate is %e, ep5=%e,sal=%e,sf=%e,sr1=%e,sr2=%e\n",
//                     meltrate, ep5, sal, sf, sr1, sr2); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POMeltingParam3eqn::adlprt(PetscReal salz,PetscReal temp_insitu, PetscReal pres, PetscReal &adlprt_out){
// Berechnet aus dem Salzgehalt/psu (SALZ), der in-situ Temperatur/degC
// (TEMP) und dem in-situ Druck/dbar (PRES) den adiabatischen Temperatur-
// gradienten/(K Dbar^-1) ADLPRT.
// Checkwert: ADLPRT =     3.255976E-4 K dbar^-1
//       fuer SALZ   =    40.0 psu
//            TEMP   =    40.0 DegC
//            PRES   = 10000.000 dbar

  PetscReal ds;
  const PetscReal s0 = 35.0;
  const PetscReal a0 = 3.5803e-5,   a1 = 8.5258e-6,   a2 = -6.8360e-8, a3 = 6.6228e-10;
  const PetscReal b0 = 1.8932e-6,   b1 = -4.2393e-8;
  const PetscReal c0 = 1.8741e-8,   c1 = -6.7795e-10, c2 = 8.7330e-12, c3 = -5.4481e-14;
  const PetscReal d0 = -1.1351e-10, d1 = 2.7759e-12;
  const PetscReal e0 = -4.6206e-13, e1 = 1.8676e-14,  e2 = -2.1687e-16;

  ds = salz-s0;
  adlprt_out = (( ( (e2*temp_insitu + e1)*temp_insitu + e0 )*pres + ( (d1*temp_insitu + d0)*ds
               + ( (c3*temp_insitu + c2)*temp_insitu + c1 )*temp_insitu + c0 ) )*pres
               + (b1*temp_insitu + b0)*ds +  ( (a3*temp_insitu + a2)*temp_insitu + a1 )*temp_insitu + a0);
  
  return 0;
}

PetscErrorCode POMeltingParam3eqn::pttmpr(PetscReal salz,PetscReal temp_insitu,PetscReal pres,PetscReal rfpres,
                                          PetscReal& thetao){
// Berechnet aus dem Salzgehalt/psu (SALZ), der in-situ Temperatur/degC
// (TEMP) und dem in-situ Druck/dbar (PRES) die potentielle Temperatur/
// degC (PTTMPR) bezogen auf den Referenzdruck/dbar (RFPRES). Es wird
// ein Runge-Kutta Verfahren vierter Ordnung verwendet.
// Checkwert: PTTMPR = 36.89073 DegC
//       fuer SALZ   =    40.0 psu
//            TEMP   =    40.0 DegC
//            PRES   = 10000.000 dbar
//            RFPRES =     0.000 dbar

  PetscReal ct2  = 0.29289322 , ct3  = 1.707106781;
  PetscReal cq2a = 0.58578644 , cq2b = 0.121320344;
  PetscReal cq3a = 3.414213562, cq3b = -4.121320344;

  //real salz,temp_insitu,pres,rfpres
  PetscReal p,t,dp,dt,q, dd;
  //real adlprt

  p  = pres;
  t  = temp_insitu;
  dp = rfpres-pres;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  t  = t +0.5*dt;
  q = dt;
  p  = p +0.5*dp;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  t  = t + ct2*(dt-q);
  q  = cq2a*dt + cq2b*q;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  t  = t + ct3*(dt-q);
  q  = cq3a*dt + cq3b*q;
  p  = rfpres;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  thetao = t+ (dt-q-q)/6.0;
  
  return 0;
}

PetscErrorCode POMeltingParam3eqn::potit(PetscReal salz,PetscReal thetao,PetscReal pres,PetscReal rfpres, PetscReal &temp_insitu_out){
  // *********************************************************************
  // Berechnet aus dem Salzgehalt[psu] (SALZ), der pot. Temperatur[oC]
  // (PT) und dem Referenzdruck[dbar] (REFPRES) die in-situ Temperatur
  // [oC] (TIN) bezogen auf den in-situ Druck[dbar] (PRES) mit Hilfe
  // eines Iterationsverfahrens aus.

  //integer iter
  //real salz,thetao,pres,rfpres,tin
  //real epsi,tpmd,pt1,ptd,pttmpr

  //tpmd / 0.001 /
  PetscReal tpmd = 0.001, epsi = 0., tin, pt1, ptd;
  
  for (PetscInt iter=0; iter<101; ++iter){
  //do iter=1,100
    tin  = thetao+epsi;
//     print "tin=" + str(tin)
    pttmpr(salz,tin,pres,rfpres,pt1);
//     PetscErrorCode ierr = verbPrintf(2, grid.com, "pt1=%e\n", pt1); CHKERRQ(ierr);
//     print "pt1=" + str(pt1)
    ptd  = pt1-thetao;
//     print "ptd=" + str(ptd)
    if(PetscAbs(ptd) < tpmd){
      break;
    }else{
      epsi = epsi-ptd;
    }
    if(iter==100){ SETERRQ(1, "in situ temperature calculation not converging."); }
  }
  
  temp_insitu_out = tin;
  return 0;
}
