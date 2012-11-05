// Copyright (C) 2010, 2011, 2012 Constantine Khroulev
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

#ifndef __PISMProf_hh
#define __PISMProf_hh

#ifndef PISM_PROFILE
#define PISM_PROFILE
#endif

#include <string>
#include <vector>
#include <petsc.h>
#include "PISMNCFile.hh"

/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

//! \brief A class storing and writing PISM profiling event data.
class PISMEvent {
public:
  PISMEvent();
  string name,			//!< NetCDF variable name
    description,		//!< NetCDF variable long_name attribute
    units;                      //!< NetCDF variable units
  int parent;			//!< index of the parent event
  PetscLogDouble start_time;	//!< event start time
  double total_time;		//!< total time spent in an event; includes
				//!< time spent doing nothing
  PetscLogEvent petsc_event;
};

//! PISM profiler class.
/*!
  Usage example:

  \code
  PISMProf *prof;
  prof = new PISMProf(grid, config);
  int event;

  event = prof->create("event_varname", "event_description");

  prof->begin(event);
  // do stuff
  prof->end(event);

  ierr = prof->save_report("prof.nc"); CHKERRQ(ierr); 

  delete prof;
  \endcode
 */
class PISMProf {
public:
  PISMProf(MPI_Comm c, PetscMPIInt r, PetscMPIInt s);
  ~PISMProf() {}
  int create(string name, string description);
  int get(string name);
  void begin(int index);
  void end(int index);
  PetscErrorCode barrier();
  PetscErrorCode save_report(string filename);
  void set_grid_size(int n);
  int Nx, Ny;
protected:
  vector<PISMEvent> events;
  int current_event;
  PetscMPIInt rank, size;
  MPI_Comm com;

  PetscErrorCode save_report(int index, const PISMNCFile &nc, string name);
  PetscErrorCode find_variables(PISMNCFile &nc, string name, bool &exists);
  PetscErrorCode define_variable(const PISMNCFile &nc, string name);
  PetscErrorCode create_dimensions(const PISMNCFile &nc);
};

/**
   The following functions were added to provide easy access to PETSc's
   profiling features, while bypassing some of the odd PISM profiling code.
   For example, the PISM profiling stuff above may segfault if the events do not
   follow a hierarchical order (start-a, start-b, end-b, end-a is expected; 
   start-a, start-b, end-a, end-b will cause problems).  The following will
   allow this event pattern, though still will have trouble with events that 
   intersect themselves (start-a, start-a, end-a, end-a).  This is a limitation
   of PETSc's profiling, and we can live with it.

   These also get around the problem of having to pass a single instance of
   PISMProf to every bit of code that wishes to do some profiling.  

   Note: this functionality isn't currently thread safe, though it wouldn't take
   but a mutex or two to fix that.
*/

/**
   The following is a list of known event identifiers.  This is certainly not an 
   exhaustive listing, as events not listed here may be used.
*/
#define PISM_IO_EVENT "IO"
#define PISM_IO_DATASET_EVENT "IO-Dataset"
#define PISM_IO_METADATA_EVENT "IO-Metadata"
#define PISM_IO_READ_EVENT "IO-Read"
#define PISM_IO_WRITE_EVENT "IO-Write"

#define PISM_INIT_STAGE "Initialization"
#define PISM_STEPPING_STAGE "Stepping"
#define PISM_PRIMARY_OUTPUT_STAGE "PrimaryOutput"
#define PISM_SHUTDOWN_STAGE "Shutdown"

//! Starts recording the named event.
/*!
  Starts recording the named event.  If this is the first use of the event, it
  will be registered with PETSc.

  \param name A name that uniquely identifies the log event.  
  \return 0 if the operation was successful, and an error code otherwise.
 */
PetscErrorCode PISMLogEventBegin(string name);

//! Stops recording the named event.
/*!
  Stops recording the named event.  It is assumed that this event was previously
  started with a call to PISMLogEventBegin.

  \param name A name that uniquely identifies the log event.
  \return 0 if the operation was successful, and an error code otherwise.
 */
PetscErrorCode PISMLogEventEnd(string name);

//! Indicates the beginning of an application stage.
/*!
  Indicates the beginning of an application stage (like initializatin, shutdown,
  stepping, etc.).  This is basically a wrapper around the PetscLogStagePush
  function, and it will take care of the stage registration the first time this
  function is called with a given stage name.

  \param name The name of the stage that will be activated.
  \return 0 if the operation was successful, and an error code otherwise.
*/
PetscErrorCode PISMLogStagePush(string name);

//! Indicates the end of the currently active log stage.
/*!
  Indicates the end of the currently active log stage.  This is basically a
  wrapper around PetscLogStagePop().
  
  \return 0 if the operation was successful, and an error code otherwise.
*/
PetscErrorCode PISMLogStagePop();

#endif // __PISMProf_hh
