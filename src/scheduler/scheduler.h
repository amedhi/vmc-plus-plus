/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:48:00
*----------------------------------------------------------------------------*/
// File: scheduler.h 
// Definition of the Scheduler class.

#ifndef SCHEDULER_SCHEDULER_H
#define SCHEDULER_SCHEDULER_H

#include <iostream>
#include "task.h"
#include "cmdargs.h"
#include "inputparams.h"
#include <chrono>

namespace scheduler {

int start(int argc, const char *argv[], const AbstractTask& theTask);

class Scheduler 
{
public:
  //Scheduler(): simmaster(0) {};
  //Scheduler(Task& theTask) {}
  Scheduler() {}
  ~Scheduler() {}
  virtual int run(void);

protected:
  AbstractWorker* theWorker;
  input::Parameters parms;
  bool valid_{false};

private:
  //unsigned simmaster;
  //TaskParams parms;
};

class MasterScheduler : public Scheduler
{
public:
  MasterScheduler(int argc, const char *argv[], const AbstractTask& theTask);
  //MasterScheduler(int argc, const char *argv[], const Task&);
  MasterScheduler() = delete;
  ~MasterScheduler() {};
  int run(void) override;

private:
  CommandArg cmdarg;
  input::JobInput input;
  unsigned int task_size;

  std::string elapsed_time(const std::chrono::steady_clock::time_point& start_time,
    const std::chrono::steady_clock::time_point& end_time) const;
};

} // end namespace input

#endif
