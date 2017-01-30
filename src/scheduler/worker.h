/*---------------------------------------------------------------------------
* task.h: Class representing a task to be run.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 10:57:03
* Last Modified by:   amedhi
* Last Modified time: 2017-01-29 22:10:51
*----------------------------------------------------------------------------*/
#include <iostream>
#include "inputparams.h"

#ifndef SCHEDULER_WORKER_H
#define SCHEDULER_WORKER_H

namespace scheduler {

class AbstractWorker
{
public:
  virtual ~AbstractWorker() {};
  virtual int start(input::Parameters& p) = 0; // start all runs
  virtual void run(void) = 0; // run for some time (in seconds)
  virtual void finish(void) = 0; // mark as finished
  virtual void halt(void) = 0; // halt all runs, simulation is finished        
};

class Worker : public AbstractWorker
{
public:
  Worker();
  virtual ~Worker();
  
  //virtual void construct(); // needs to be called to finish construction
  int start(input::Parameters& p) override; // start simulation
  void run(void) override;// run a few steps and return control
  virtual void dostep(void)=0; // do a step
  // bool started() const { return started_;}
  void halt(void) override;
protected:
  bool finished_;

private:
  bool started_; // is the task running?
};


} // end namespace input

#endif
