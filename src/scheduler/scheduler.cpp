/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 12:44:04
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-08 21:55:35
*----------------------------------------------------------------------------*/
#include <iostream>
#include "scheduler.h"

namespace scheduler {

int start(int argc, const char *argv[], const AbstractTask& theTask)
{
  //Scheduler* theScheduler;
  //theScheduler = new MasterScheduler(argc, argv, theTask);
  //int res = theScheduler->run();
  MasterScheduler master_scheduler(argc, argv, theTask);
  int res = master_scheduler.run();
	return res;	
}

int Scheduler::run(void) 
{
  // Task task;
  // valid = input.read_params(0);
  // task.init_task_param()
  return 0;
}

} // end namespace scheduler
