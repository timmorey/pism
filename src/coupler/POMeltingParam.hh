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

#ifndef _POMELTINGPARAM_H_
#define _POMELTINGPARAM_H_

#include "PASDirectForcing.hh"
#include "PISMOcean.hh"

class POMeltingParam3eqn : public PDirectForcing<PISMOceanModel>
{
public:
  POMeltingParam3eqn(IceGrid &g, const NCConfigVariable &conf)
    : PDirectForcing<PISMOceanModel>(g, conf)
  {
    temp_name       = "oceantemp";
    mass_flux_name  = "salinity";
    bc_option_name = "-ocean_bc_file";
  }

  virtual ~POMeltingParam3eqn() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result) {
    result = sea_level;
    return 0;
  }

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode cavity_heat_water_fluxes_3eq(PetscReal temp,
               PetscReal sal, PetscReal tin, PetscReal zice, PetscReal rhow,
               PetscReal rhoi, PetscReal rho, PetscReal &meltrate);
  virtual PetscErrorCode adlprt(PetscReal salz,PetscReal temp,PetscReal pres);
  virtual PetscErrorCode pttmpr(PetscReal salz,PetscReal temp,PetscReal pres,PetscReal rfpres);
  virtual PetscErrorCode potit(PetscReal salz,PetscReal pt,PetscReal pres,PetscReal rfpres);
  protected:
    IceModelVec2S *ice_thickness, *bed; // is not owned by this class

};


#endif /* _POMELTINGPARAM_H_ */
