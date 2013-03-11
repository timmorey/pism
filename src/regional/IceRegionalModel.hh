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

#ifndef ICEREGIONALMODEL_HH_
#define ICEREGIONALMODEL_HH_


#include "iceModel.hh"

#include <petsc.h>
#include <string>


class CoarseGrid;

//! \brief A version of the PISM core class (IceModel) which knows about the
//! no_model_mask and its semantics.
class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid &g, NCConfigVariable &c, NCConfigVariable &o);

public:
  virtual PetscErrorCode attach_coarse_grid(const std::string& filename);
  virtual PetscErrorCode step(bool do_mass_continuity, bool do_energy, bool do_age, bool do_skip);

protected:
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode bootstrap_2d(string filename);
  virtual PetscErrorCode initFromFile(string filename);
  virtual PetscErrorCode model_state_setup();
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode allocate_basal_yield_stress();
  virtual PetscErrorCode massContExplicitStep();
  virtual void cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     planeStar<PISMVector2> input_velocity,
                                     planeStar<PetscScalar> input_flux,
                                     planeStar<PetscScalar> &output_velocity,
                                     planeStar<PetscScalar> &output_flux);
  virtual PetscErrorCode enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                                                 PetscScalar* bulgeCount);
  virtual PetscErrorCode setFromOptions();

private:
  IceModelVec2Int no_model_mask;
  IceModelVec2S   usurfstore, thkstore;
  IceModelVec2S   bmr_stored;
  PetscErrorCode  set_no_model_strip(PetscReal stripwidth);

  CoarseGrid* coarse_grid;
};


#endif /* ICEREGIONALMODEL_HH_ */
