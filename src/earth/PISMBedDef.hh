// Copyright (C) 2010 Constantine Khroulev
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

#ifndef __PISMBedDef_hh
#define __PISMBedDef_hh

#include "../base/PISMComponent.hh"
#include "../base/iceModelVec.hh"
#include "deformation.hh"

//! PISM bed deformation model (base class).
/*! Unlike other PISMComponent derived classes, the update() method of
  PISMBedDef has side-effects (modifies IceModel data memebers).
 */
class PISMBedDef : public PISMComponent {
public:
  PISMBedDef(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMBedDef() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years) = 0;
protected:
  PetscErrorCode pismbeddef_allocate(); // packaged to simplify error checking
  PetscErrorCode compute_uplift(PetscScalar dt_beddef);
  PetscReal t_beddef_last;		//!< last bed deformation update year

  IceModelVec2S topg_last;	//!< last bed elevation
  IceModelVec2S *thk,		//!< pointer to the current ice thickness
    *topg,			//!< pointer to the current bed elevation
    *uplift;			//!< pointer to the bed uplift rate field
};

//! Pointwide isostasy bed deformation model.
class PBPointwiseIsostasy : public PISMBedDef {
public:
  PBPointwiseIsostasy(IceGrid &g, const NCConfigVariable &conf); 
  virtual ~PBPointwiseIsostasy() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
protected:
  PetscErrorCode allocate();
  IceModelVec2S thk_last;	//!< last ice thickness
};

#if (PISM_HAVE_FFTW==1)
#include <fftw3.h>

//! A wrapper class around BedDeformLC.
class PBLingleClark : public PISMBedDef {
public:
  PBLingleClark(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PBLingleClark();

  PetscErrorCode init(PISMVars &vars);
  PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
protected:
  PetscErrorCode allocate();
  PetscErrorCode deallocate();
  PetscErrorCode transfer_to_proc0(IceModelVec2S *source, Vec result);
  PetscErrorCode transfer_from_proc0(Vec source, IceModelVec2S *result);
  Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
  VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
  // Vecs on processor 0:
  Vec Hp0,			//!< ice thickness
    bedp0,			//!< bed elevation
    Hstartp0,			//!< initial (start-of-the-run) thickness
    bedstartp0,			//!< initial bed elevation
    upliftp0;			//!< bed uplift
  BedDeformLC bdLC;
};
#endif	// PISM_HAVE_FFTW

#endif	// __PISMBedDef_hh