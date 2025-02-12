// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "rerun_min_drude.h"

#include "domain.h"
#include "error.h"
#include "finish.h"
#include "integrate.h"
#include "modify.h"
#include "output.h"
#include "read_dump.h"
#include "timer.h"
#include "update.h"
#include "min.h"
#include "comm.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RerunMinDrude::RerunMinDrude(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void RerunMinDrude::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Rerun command before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal rerun command");

  // list of dump files = args until a keyword

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) break;
    if (strcmp(arg[iarg],"last") == 0) break;
    if (strcmp(arg[iarg],"every") == 0) break;
    if (strcmp(arg[iarg],"skip") == 0) break;
    if (strcmp(arg[iarg],"start") == 0) break;
    if (strcmp(arg[iarg],"stop") == 0) break;
    if (strcmp(arg[iarg],"dump") == 0) break;
    if (strcmp(arg[iarg],"post") == 0) break;
    iarg++;
  }
  int nfile = iarg;
  if (nfile == 0 || nfile == narg) error->all(FLERR,"Illegal rerun command");

  // parse optional args up until "dump"
  // use MAXBIGINT -1 so Output can add 1 to it and still be a big int

  bigint first = 0;
  bigint last = MAXBIGINT - 1;
  int nevery = 0;
  int nskip = 1;
  int startflag = 0;
  int stopflag = 0;
  int postflag = 0;
  bigint start = -1;
  bigint stop = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      first = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (first < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"last") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      last = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (last < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"skip") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nskip = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nskip <= 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      startflag = 1;
      start = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (start < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"stop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      stopflag = 1;
      stop = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (stop < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      postflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dump") == 0) {
      break;
    } else error->all(FLERR,"Illegal rerun command");
  }

  int nremain = narg - iarg - 1;
  if (nremain <= 0) error->all(FLERR,"Illegal rerun command");
  if (first > last) error->all(FLERR,"Illegal rerun command");
  if (startflag && stopflag && start > stop)
    error->all(FLERR,"Illegal rerun command");

  // pass list of filenames to ReadDump
  // along with post-"dump" args and post-"format" args

  auto rd = new ReadDump(lmp);

  rd->store_files(nfile,arg);
  if (nremain)
    nremain = rd->fields_and_keywords(nremain,&arg[narg-nremain]);
  else nremain = rd->fields_and_keywords(0,nullptr);
  if (nremain) rd->setup_reader(nremain,&arg[narg-nremain]);
  else rd->setup_reader(0,nullptr);

  // perform the pseudo run
  // invoke lmp->init() only once
  // read all relevant snapshots
  // use setup_minimal() since atoms are already owned by correct procs
  // addstep_compute_all() insures energy/virial computed on every snapshot

  update->whichflag = 2;

  if (startflag) update->beginstep = update->firststep = start;
  else update->beginstep = update->firststep = first;
  if (stopflag) update->endstep = update->laststep = stop;
  else update->endstep = update->laststep = last;

  int firstflag = 1;
  int ndump = 0;

  lmp->init();
  bigint ntimestep = rd->seek(first,0);
  if (ntimestep < 0)
    error->all(FLERR,"Rerun dump file does not contain requested snapshot");

  int count = 0;
  while (true) {
    ndump++;
    rd->header(firstflag);
    update->reset_timestep(ntimestep, false);
    rd->atoms();

    update->etol = 1e-10;
    update->ftol = 1e-10;
    update->nsteps = 1000;
    update->max_eval = 10000;
    
    int current_step = update->ntimestep;
    update->vflag_global =  current_step;

    update->beginstep = update->firststep = update->ntimestep;
    update->endstep = update->firststep + update->nsteps;
  
    timer->init_timeout();
    if (count == 0){
      update->minimize->setup_drude();
    } else {
      update->minimize->setup_minimal(1);
    } 

    timer->init();
    timer->barrier_start();
    update->minimize->run(update->nsteps);
    timer->barrier_stop();

    output->next_dump_any = ntimestep;
    if (firstflag) output->setup();
    else if (output->next) output->write(ntimestep);

    update->minimize->force_clear();
    update->minimize->niter = 0;
    update->minimize->neval = 0;

    firstflag = 0;
    ntimestep = rd->next(ntimestep,last,nevery,nskip);
    if (stopflag && ntimestep > stop)
      error->all(FLERR,"Read rerun dump file timestep > specified stop");
    if (ntimestep < 0) break;

    count++;
  }

  // insure thermo output on last dump timestep

  output->next_thermo = update->ntimestep;
  output->write(update->ntimestep);

  timer->barrier_stop();

  update->minimize->cleanup();

  // set update->nsteps to ndump for Finish stats to print

  update->nsteps = ndump;

  Finish finish(lmp);
  finish.end(postflag);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;

  // clean-up

  delete rd;
}
