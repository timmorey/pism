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

#include "PISMProf.hh"
#include "PISMNC3File.hh"
#include "pism_const.hh"

#include <map>

/// PISMEvent
PISMEvent::PISMEvent() {
  name = "unknown";
  units = "seconds";
  start_time = -1;
  total_time = 0;
  parent = -1;
  petsc_event = 0;
}

/// PISMProf
PISMProf::PISMProf(MPI_Comm c, PetscMPIInt r, PetscMPIInt s) {
  com = c;
  rank = r;
  size = s;
  Nx = s;
  Ny = 1;
  current_event = -1;

  PISMEvent tmp;
  tmp.name = "processor_rank";
  tmp.description = "processor rank";
  tmp.total_time = r;
  tmp.units = "count";
  events.push_back(tmp);
}

void PISMProf::set_grid_size(int n) {

  PISMEvent tmp;
  tmp.name = "processor_grid_size";
  tmp.description = "number of map-plane grid points in a processor's subdomain";
  tmp.total_time = n;
  tmp.units = "count";
  events.push_back(tmp);
}


# define PismLogEventRegister(name,cookie,event) PetscLogEventRegister((name),(cookie),(event))

//! Create a profiling event.
/*!
 * Checks if an event with this name already exists.
 */
int PISMProf::create(string name, string description) {
  PISMEvent tmp;
  int index = get(name);

  if (index != -1)
    return index;

  tmp.name = name;
  tmp.description = description;
  PismLogEventRegister(name.c_str(), 0, &tmp.petsc_event);

  events.push_back(tmp);

  return (int)events.size() - 1;
}

//! \brief Get an integer (index) corresponding to an event.
/*!
 * Returns -1 if an event was not found.
 */
int PISMProf::get(string name) {

  for (unsigned int i = 0; i < events.size(); ++i)
    if (events[i].name == name)
      return i;

  return -1;
}

// This ensures that begin() and end() code is optimized away if
// PISM_PROFILE is not defined.
#ifndef PISM_PROFILE
void PISMProf::begin(int) {}
#else
//! Begin a profiling event.
void PISMProf::begin(int index) {
  PISMEvent &event = events[index];

  // fprintf(stderr, "Rank %d: Start '%s'\n", rank, events[index].description.c_str());

  event.parent = current_event;

  PetscGetTime(&event.start_time);
  PetscLogEventBegin(event.petsc_event, 0, 0, 0, 0);

  current_event = index;
}
#endif

#ifndef PISM_PROFILE
void PISMProf::end(int) {}
#else
//! End a profiling event.
void PISMProf::end(int /*index*/) {
  PetscLogDouble time;
  PISMEvent &event = events[current_event];

  PetscGetTime(&time);
  PetscLogEventEnd(event.petsc_event, 0, 0, 0, 0);

  event.total_time += (time - event.start_time);

  current_event = event.parent;
}
#endif

//! Save a profiling report to a file.
PetscErrorCode PISMProf::save_report(string filename) {
  PetscErrorCode ierr;
  PISMNC3File nc(com, rank);

  ierr = nc.create(filename); CHKERRQ(ierr);

  ierr = nc.def_dim("x", Nx); CHKERRQ(ierr);
  ierr = nc.def_dim("y", Ny); CHKERRQ(ierr);

  for (unsigned int j = 0; j < events.size(); ++j) {
    string &name = events[j].name;
    if (name == "unknown")
      continue;

    double time = events[j].total_time,
      max_time = 0;

    PISMGlobalMax(&time, &max_time, com);

    // ignore events that took less than 0.001 seconds
    if (max_time < 1e-3)
      continue;

    ierr = save_report(j, nc, name); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Save the report corresponding to an event events[index].
PetscErrorCode PISMProf::save_report(int index, const PISMNCFile &nc, string variable_name) {
  PetscErrorCode ierr;
  double data[1];
  vector<unsigned int> start(2), count(2);

  ierr = define_variable(nc, variable_name); CHKERRQ(ierr);

  int parent_index = events[index].parent;
  string descr = events[index].description,
    parent = parent_index == -1 ? "root" : events[parent_index].name;

  ierr = nc.put_att_text(variable_name, "units",     events[index].units); CHKERRQ(ierr); 
  ierr = nc.put_att_text(variable_name, "parent",    parent); CHKERRQ(ierr); 
  ierr = nc.put_att_text(variable_name, "long_name", descr); CHKERRQ(ierr);

  data[0] = events[index].total_time;

  start[0] = rank % (Ny);
  start[1] = rank / (Ny);

  count[0] = 1;  count[1] = 1;

  ierr = nc.enddef(); CHKERRQ(ierr);
  ierr = nc.put_vara_double(variable_name, start, count, data); CHKERRQ(ierr);

  return 0;
}

//! Define a NetCDF variable to store a profiling report in.
PetscErrorCode PISMProf::define_variable(const PISMNCFile &nc, string name) {
  PetscErrorCode ierr;
  vector<string> dims;
  bool exists = false;

  ierr = nc.inq_varid(name, exists); CHKERRQ(ierr);

  if (exists)
    return 0;

  dims.push_back("y"); dims.push_back("x");

  ierr = nc.redef(); CHKERRQ(ierr);
  ierr = nc.def_var(name, PISM_DOUBLE, dims); CHKERRQ(ierr);
  ierr = nc.enddef(); CHKERRQ(ierr);

  return 0;
}

//! Barrier synchroniation. (Calls MPI_Barrier.)
PetscErrorCode PISMProf::barrier() {

#ifdef PISM_PROFILE
  PetscErrorCode ierr;
  ierr = MPI_Barrier(com); CHKERRQ(ierr);
#endif

  return 0;
}

/**
   The following is the implementation of the profiling functions defined at the
   bottom of PISMProf.hh.  We use global variables here to ensure that we have
   exactly one instance of these per process.
 */

PetscClassId g_DefaultEventClassId = -1;
std::map<string, PetscLogEvent> g_RegisteredEvents;
std::map<string, PetscLogStage> g_RegisteredStages;

PetscErrorCode PISMLogEventBegin(string name) {
  PetscErrorCode retval = 0;
  std::map<string, PetscLogEvent>::iterator iter;
  PetscLogEvent event = -1;

  iter = g_RegisteredEvents.find(name);
  if(iter == g_RegisteredEvents.end()) {
    // The event has not yet been registered, so we must do that first.

    if(0 > g_DefaultEventClassId) {
      // We have not even registered our class id yet, and we must do that
      // before we register the event.
      retval = PetscClassIdRegister("PISM-DefaultClassId", &g_DefaultEventClassId); 
      CHKERRQ(retval);
    }

    PetscLogEventRegister(name.c_str(), g_DefaultEventClassId, &event); 
    CHKERRQ(retval);

    g_RegisteredEvents[name] = event;
  
  } else {
    event = iter->second;
  }

  if(0 > event) {
    // Something went wrong, and we failed to register the event.
    fprintf(stderr, "Failed to register event \"%s\"\n", name.c_str());
    retval = -1;
  } else {
    retval = PetscLogEventBegin(event, 0, 0, 0, 0);  CHKERRQ(retval);
  }

  return retval;
}

PetscErrorCode PISMLogEventEnd(string name) {
  PetscErrorCode retval = 0;
  std::map<string, PetscLogEvent>::iterator iter;

  iter = g_RegisteredEvents.find(name);
  if(iter != g_RegisteredEvents.end()) {
    retval = PetscLogEventEnd(iter->second, 0, 0, 0, 0);  CHKERRQ(retval);
  } else {
    // The event is not registered, so we can't end it.
    fprintf(stderr, "Event \"%s\" is not registered.\n", name.c_str());
    retval = -1;
  }

  return retval;
}

PetscErrorCode PISMLogStagePush(string name) {
  PetscErrorCode retval = 0;
  std::map<string, PetscLogStage>::iterator iter;
  PetscLogStage stage = -1;

  iter = g_RegisteredStages.find(name);
  if(iter == g_RegisteredStages.end()) {
    PetscLogStageRegister(name.c_str(), &stage);  CHKERRQ(retval);
    g_RegisteredStages[name] = stage;
  } else {
    stage = iter->second;
  }

  if(0 > stage) {
    // Something went wrong, and we failed to register the stage.
    fprintf(stderr, "Failed to register stage \"%s\"\n", name.c_str());
    retval = -1;
  } else {
    retval = PetscLogStagePush(stage);  CHKERRQ(retval);
  }

  return retval;
}

PetscErrorCode PISMLogStagePop() {
  PetscErrorCode retval = 0;
  retval = PetscLogStagePop();  CHKERRQ(retval);
  return retval;
}
