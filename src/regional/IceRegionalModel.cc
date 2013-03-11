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


#include "CoarseGrid.hh"
#include "IceGrid.hh"
#include "IceRegionalModel.hh"
#include "pism_options.hh"
#include "PISMConstantYieldStress.hh"
#include "regional.hh"


IceRegionalModel::IceRegionalModel(IceGrid &g, NCConfigVariable &c, NCConfigVariable &o)
  : IceModel(g,c,o),
    coarse_grid(0) {

}

PetscErrorCode IceRegionalModel::attach_coarse_grid(const std::string& filename) {

  PetscErrorCode retval = 0;

  if(coarse_grid) {
    delete coarse_grid;
    coarse_grid = 0;
  }

  coarse_grid = new CoarseGrid(filename);
  retval = coarse_grid->SetAreaOfInterest(grid.x[grid.xs], grid.x[grid.xs + grid.xm],
                                          grid.y[grid.ys], grid.y[grid.ys + grid.ym]);
  CHKERRQ(retval);

  return retval;
}

PetscErrorCode IceRegionalModel::step(bool do_mass_continuity, bool do_energy, bool do_age, bool do_skip) {
  PetscErrorCode retval = 0;
  double value;

  if(this->grid.time->current() >= this->grid.time->end()) {
    // TODO: Check to see if this will be the last time step.  If so, then we may
    // need to write out the current (next to last) time step so that we can
    // interpolate the model state into a coarse model.
  }

  if(coarse_grid) {
    this->no_model_mask.begin_access();
    this->vH.begin_access();

    PetscReal t = this->grid.time->current();
    for (PetscInt i = grid.xs; i < grid.xs+grid.xm; ++i) {
      PetscReal x = this->grid.x[i];
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if (no_model_mask(i, j) > 0.5) {
          PetscReal y = this->grid.y[i];

          this->coarse_grid->Interpolate(this->vH.string_attr("name", 0), x, y, 0.0, t, &value);
          this->vH(i, j) = value;
        }
      }
    }

    this->no_model_mask.end_access();
    this->vH.end_access();
  }

  retval = IceModel::step(do_mass_continuity, do_energy, do_age, do_skip);
  CHKERRQ(retval);

  return retval;
}

//! \brief
PetscErrorCode IceRegionalModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);

  ierr = config.flag_from_option("ssa_dirichlet_bc", "ssa_dirichlet_bc"); CHKERRQ(ierr);

  return 0;
}

//! \brief Set no_model_mask variable to have value 1 in strip of width 'strip'
//! m around edge of computational domain, and value 0 otherwise.
PetscErrorCode IceRegionalModel::set_no_model_strip(PetscReal strip) {
  PetscErrorCode ierr;

    ierr = no_model_mask.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if (grid.x[i] <= grid.x[0]+strip || grid.x[i] >= grid.x[grid.Mx-1]-strip) {
          no_model_mask(i, j) = 1;
        } else if (grid.y[j] <= grid.y[0]+strip || grid.y[j] >= grid.y[grid.My-1]-strip) {
          no_model_mask(i, j) = 1;
        } else {
          no_model_mask(i, j) = 0;
        }
      }
    }
    ierr = no_model_mask.end_access(); CHKERRQ(ierr);

    ierr = no_model_mask.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);

    ierr = no_model_mask.update_ghosts(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "  creating IceRegionalModel vecs ...\n"); CHKERRQ(ierr);

  // stencil width of 2 needed for surfaceGradientSIA() action
  ierr = no_model_mask.create(grid, "no_model_mask", true, 2); CHKERRQ(ierr);
  ierr = no_model_mask.set_attrs("model_state", // ensures that it gets written at the end of the run
    "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)",
    "", ""); CHKERRQ(ierr); // no units and no standard name
  double NMMASK_NORMAL   = 0.0,
         NMMASK_ZERO_OUT = 1.0;
  vector<double> mask_values(2);
  mask_values[0] = NMMASK_NORMAL;
  mask_values[1] = NMMASK_ZERO_OUT;
  ierr = no_model_mask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = no_model_mask.set_attr("flag_meanings", "normal special_treatment"); CHKERRQ(ierr);
  no_model_mask.output_data_type = PISM_BYTE;
  no_model_mask.time_independent = true;
  ierr = no_model_mask.set(NMMASK_NORMAL); CHKERRQ(ierr);
  ierr = variables.add(no_model_mask); CHKERRQ(ierr);

  // stencil width of 2 needed for differentiation because GHOSTS=1
  ierr = usurfstore.create(grid, "usurfstore", true, 2); CHKERRQ(ierr);
  ierr = usurfstore.set_attrs(
    "model_state", // ensures that it gets written at the end of the run
    "saved surface elevation for use to keep surface gradient constant in no_model strip",
    "m",
    ""); CHKERRQ(ierr); //  no standard name
  ierr = variables.add(usurfstore); CHKERRQ(ierr);

  // stencil width of 1 needed for differentiation
  ierr = thkstore.create(grid, "thkstore", true, 1); CHKERRQ(ierr);
  ierr = thkstore.set_attrs(
    "model_state", // ensures that it gets written at the end of the run
    "saved ice thickness for use to keep driving stress constant in no_model strip",
    "m",
    ""); CHKERRQ(ierr); //  no standard name
  ierr = variables.add(thkstore); CHKERRQ(ierr);

  // Note that the name of this variable (bmr_stored) does not matter: it is
  // *never* read or written. We make a copy of bmelt instead.
  ierr = bmr_stored.create(grid, "bmr_stored", true, 2); CHKERRQ(ierr);
  ierr = bmr_stored.set_attrs("internal",
                              "time-independent basal melt rate in the no-model-strip",
                              "m s-1", ""); CHKERRQ(ierr);

  if (config.get_flag("ssa_dirichlet_bc")) {
    // remove the bcflag variable from the dictionary
    variables.remove("bcflag");

    ierr = variables.add(no_model_mask, "bcflag"); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceRegionalModel::model_state_setup() {
  PetscErrorCode ierr;

  ierr = IceModel::model_state_setup(); CHKERRQ(ierr);

  // Now save the basal melt rate at the beginning of the run.
  ierr = bmr_stored.copy_from(vbmr); CHKERRQ(ierr);

  bool zgwnm;
  ierr = PISMOptionsIsSet("-zero_grad_where_no_model", zgwnm); CHKERRQ(ierr);
  if (zgwnm) {
    ierr = thkstore.set(0.0); CHKERRQ(ierr);
    ierr = usurfstore.set(0.0); CHKERRQ(ierr);
  }

  bool nmstripSet;
  PetscReal stripkm = 0.0;
  ierr = PISMOptionsReal("-no_model_strip",
                         "width in km of strip near boundary in which modeling is turned off",
       stripkm, nmstripSet);

  if (nmstripSet) {
    ierr = verbPrintf(2, grid.com,
                      "* Option -no_model_strip read... setting boundary strip width to %.2f km\n",
                      stripkm); CHKERRQ(ierr);
    ierr = set_no_model_strip(convert(stripkm, "km", "m")); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceRegionalModel::allocate_stressbalance() {
  PetscErrorCode ierr;

  bool use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_sia = config.get_flag("do_sia");

  ShallowStressBalance *my_stress_balance;
  SSB_Modifier *modifier;
  if (do_sia) {
    modifier = new SIAFD_Regional(grid, *EC, config);
  } else {
    modifier = new SSBM_Trivial(grid, *EC, config);
  }

  if (use_ssa_velocity) {
    my_stress_balance = new SSAFD_Regional(grid, *basal, *EC, config);
  } else {
    my_stress_balance = new SSB_Trivial(grid, *basal, *EC, config);
  }

  // ~PISMStressBalance() will de-allocate my_stress_balance and modifier.
  stress_balance = new PISMStressBalance(grid, my_stress_balance,
                                         modifier, ocean, config);

  // Note that in PISM stress balance computations are diagnostic, i.e. do not
  // have a state that changes in time. This means that this call can be here
  // and not in model_state_setup() and we don't need to re-initialize after
  // the "diagnostic time step".
  ierr = stress_balance->init(variables); CHKERRQ(ierr);

  if (config.get_flag("include_bmr_in_continuity")) {
    ierr = stress_balance->set_basal_melt_rate(&vbmr); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceRegionalModel::allocate_basal_yield_stress() {
  PetscErrorCode ierr;

  if (basal_yield_stress != NULL)
    return 0;

  bool use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_blatter = config.get_flag("do_blatter");

  if (use_ssa_velocity || do_blatter) {
    bool hold_tauc;
    ierr = PISMOptionsIsSet("-hold_tauc", hold_tauc); CHKERRQ(ierr);

    if (hold_tauc) {
      basal_yield_stress = new PISMConstantYieldStress(grid, config);
    } else {
      basal_yield_stress = new PISMRegionalDefaultYieldStress(grid, config, subglacial_hydrology);
    }
  }

  return 0;
}


PetscErrorCode IceRegionalModel::bootstrap_2d(string filename) {
  PetscErrorCode ierr;

  ierr = IceModel::bootstrap_2d(filename); CHKERRQ(ierr);

  ierr = usurfstore.regrid(filename, 0.0); CHKERRQ(ierr);
  ierr =   thkstore.regrid(filename, 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceRegionalModel::initFromFile(string filename) {
  PetscErrorCode  ierr;
  PIO nc(grid, "guess_mode");

  bool no_model_strip_set;
  ierr = PISMOptionsIsSet("-no_model_strip", "No-model strip, in km",
                          no_model_strip_set); CHKERRQ(ierr);

  if (no_model_strip_set) {
    ierr = no_model_mask.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com,
                    "* Initializing IceRegionalModel from NetCDF file '%s'...\n",
                    filename.c_str()); CHKERRQ(ierr);

  // Allow re-starting from a file that does not contain u_ssa_bc and v_ssa_bc.
  // The user is probably using -regrid_file to bring in SSA B.C. data.
  if (config.get_flag("ssa_dirichlet_bc")) {
    bool u_ssa_exists, v_ssa_exists;

    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_var("u_ssa_bc", u_ssa_exists); CHKERRQ(ierr);
    ierr = nc.inq_var("v_ssa_bc", v_ssa_exists); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    if (! (u_ssa_exists && v_ssa_exists)) {
      ierr = vBCvel.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
      ierr = verbPrintf(2, grid.com,
                        "PISM WARNING: u_ssa_bc and/or v_ssa_bc not found in %s. Setting them to zero.\n"
                        "              This may be overridden by the -regrid_file option.\n",
                        filename.c_str()); CHKERRQ(ierr);

      ierr = vBCvel.set(0.0); CHKERRQ(ierr);
    }
  }

  bool zgwnm;
  ierr = PISMOptionsIsSet("-zero_grad_where_no_model", zgwnm); CHKERRQ(ierr);
  if (zgwnm) {
    ierr = thkstore.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
    ierr = usurfstore.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
  }

  ierr = IceModel::initFromFile(filename); CHKERRQ(ierr);

  if (config.get_flag("ssa_dirichlet_bc")) {
      ierr = vBCvel.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  }

  if (zgwnm) {
    ierr = thkstore.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
    ierr = usurfstore.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceRegionalModel::set_vars_from_options() {
  PetscErrorCode ierr;
  bool nmstripSet;

  // base class reads the -boot_file option and does the bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

  ierr = PISMOptionsIsSet("-no_model_strip",
                          "width in km of strip near boundary in which modeling is turned off",
                          nmstripSet);
  if (!nmstripSet) {
    ierr = PetscPrintf(grid.com,
      "PISMO ERROR: option '-no_model_strip X' (X in km) is REQUIRED if '-i' is not used.\n"
      "   pismo has no well-defined semantics without it!  ENDING ...\n\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  if (config.get_flag("do_cold_ice_methods")) {
    PetscPrintf(grid.com, "PISM ERROR: pismo does not support the 'cold' mode.\n");
    PISMEnd();
  }

  return 0;
}

PetscErrorCode IceRegionalModel::massContExplicitStep() {
  PetscErrorCode ierr;

  // This ensures that no_model_mask is available in
  // IceRegionalModel::cell_interface_fluxes() below.
  ierr = no_model_mask.begin_access(); CHKERRQ(ierr);

  ierr = IceModel::massContExplicitStep(); CHKERRQ(ierr);

  ierr = no_model_mask.end_access(); CHKERRQ(ierr);

  return 0;
}

void IceRegionalModel::cell_interface_fluxes(bool dirichlet_bc,
                                             int i, int j,
                                             planeStar<PISMVector2> input_velocity,
                                             planeStar<PetscScalar> input_flux,
                                             planeStar<PetscScalar> &output_velocity,
                                             planeStar<PetscScalar> &output_flux) {

  IceModel::cell_interface_fluxes(dirichlet_bc, i, j,
                                  input_velocity,
                                  input_flux,
                                  output_velocity,
                                  output_flux);

  planeStar<int> nmm = no_model_mask.int_star(i,j);
  PISM_Direction dirs[4] = {North, East, South, West};

  for (int n = 0; n < 4; ++n) {
    PISM_Direction direction = dirs[n];

      if ((nmm.ij == 1) ||
          (nmm.ij == 0 && nmm[direction] == 1)) {
      output_velocity[direction] = 0.0;
      output_flux[direction] = 0.0;
    }
  }
  //
}

PetscErrorCode IceRegionalModel::enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                                                         PetscScalar* bulgeCount) {
  PetscErrorCode ierr;
  PetscScalar *new_enthalpy, *old_enthalpy;

  ierr = IceModel::enthalpyAndDrainageStep(vertSacrCount, liquifiedVol, bulgeCount); CHKERRQ(ierr);

  // note that the call above sets vWork3d; ghosts are comminucated later (in
  // IceModel::energyStep()).
  ierr = no_model_mask.begin_access(); CHKERRQ(ierr);

  ierr = vWork3d.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (no_model_mask(i, j) < 0.5)
        continue;

      ierr = vWork3d.getInternalColumn(i, j, &new_enthalpy); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i, j, &old_enthalpy); CHKERRQ(ierr);

      for (int k = 0; k < grid.Mz; ++k)
        new_enthalpy[k] = old_enthalpy[k];
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);

  // set vbmr; ghosts are comminucated later (in IceModel::energyStep()).
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = bmr_stored.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (no_model_mask(i, j) < 0.5)
        continue;

      vbmr(i, j) = bmr_stored(i, j);
    }
  }
  ierr = bmr_stored.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  ierr = no_model_mask.end_access(); CHKERRQ(ierr);

  return 0;
}


