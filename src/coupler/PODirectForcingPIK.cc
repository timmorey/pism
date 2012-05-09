// Copyright (C) 2011 Constantine Khroulev
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

#include "PODirectForcingPIK.hh"

PetscErrorCode PODirectForcingPIK::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean model reading ocean temperature\n"
                    "  and dummy sub-shelf mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing",
                        "absolute ocean temperature at ice shelf base",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                       "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                       "m s-1", ""); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
                    "    reading boundary conditions from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(1, "ERROR: ice thickness is not available"); }

  return 0;
}

PetscErrorCode PODirectForcingPIK::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr = update_internal(t_years, dt_years); CHKERRQ(ierr);

  if (enable_time_averaging) {
//     ierr = mass_flux.average(t, dt); CHKERRQ(ierr);
    ierr = temp.average(t, dt); CHKERRQ(ierr); 
  } else {
//     ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
    ierr = temp.get_record_years(t); CHKERRQ(ierr);
  }

  return 0;
}

// PetscErrorCode PODirectForcingPIK::shelf_base_temperature(IceModelVec2S &result) {
//   PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
//   return 0;
// }

// PetscErrorCode PODirectForcingPIK::shelf_base_mass_flux(IceModelVec2S &result) {
//   PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
//   return 0;
// }

PetscErrorCode PODirectForcingPIK::shelf_base_temperature(IceModelVec2S &result) {
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
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME task #7297
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

//! Computes mass flux in ice-equivalent m s-1, from assumption that basal heat flux rate converts to mass flux.
PetscErrorCode PODirectForcingPIK::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  PetscReal L = config.get("water_latent_heat_fusion"),
        rho_ocean = config.get("sea_water_density"),
        rho_ice = config.get("ice_density");
        //beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity");
      // K/m      Clausius-Clapeyron gradient
  const PetscScalar c_p_ocean   = 3974.0,   // J/(K*kg), specific heat capacity of ocean mixed layer
            gamma_T   = 1e-4;     // m/s, thermal exchange velocity //FIXME gamma_T should be a function of the friction velocity, not a const

  PetscScalar T_ocean_offset = 0.0;
//     T_ocean_offset = 273.15 + T_water_offset;

  // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
  // PetscReal meltrate = config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * rho_ice); // m s-1

  PetscReal meltfactor = 5e-3;
  bool meltfactorSet;
  double meltfactor_pik;
  ierr = PISMOptionsReal("-meltfactor_pik",
                           "Uses as a meltfactor as in sub-shelf-melting parameterization of martin_winkelmann11",
                           meltfactor_pik, meltfactorSet); CHKERRQ(ierr);
  if (meltfactorSet) {
    meltfactor = meltfactor_pik; //default is 5e-3 as in martin_winkelmann11
  }
//   ierr = verbPrintf(2, grid.com,"meltfactor=%f\n",meltfactor); CHKERRQ(ierr);

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = temp.begin_access(); CHKERRQ(ierr);

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        // compute T_f[i][j] according to beckmann_goosse03, which has the meaning of the freezing temperature of the ocean water directly under the shelf, (of salinity 35psu) [this is related to the Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for details]
        PetscScalar shelfbaseelev = - (rho_ice / rho_ocean) * H[i][j],
        T_f= 273.15+ (0.0939 -0.057 * 35.0 + 7.64e-4 * shelfbaseelev); // add 273.15 to get it in Kelvin... 35 is the salinity

          // compute ocean_heat_flux according to beckmann_goosse03
          // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
        PetscScalar oceanheatflux = meltfactor * rho_ocean * c_p_ocean * gamma_T * (T_ocean_offset + temp(i,j) - T_f);

          // shelfbmassflux is positive if ice is freezing on; here it is always negative:
          // same sign as OceanHeatFlux... positive if massflux FROM ice TO ocean
          //result(i,j) = oceanheatflux / (L * rho_ice) * secpera; // m a-1
          result(i,j) = oceanheatflux / (L * rho_ice); // m s-1

//         if(i==40 && j==40){
//           ierr = verbPrintf(2, grid.com,
//                   "temp=%f at i=%d, j=%d\n", temp(i,j), i,j); CHKERRQ(ierr);
//         }
      }
    }

    ierr = ice_thickness->end_access(); CHKERRQ(ierr);
    ierr = result.end_access(); CHKERRQ(ierr);
    ierr = temp.end_access(); CHKERRQ(ierr);

  return 0;
}



