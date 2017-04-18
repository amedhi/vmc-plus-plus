/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 12:44:04
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-18 15:26:55
*----------------------------------------------------------------------------*/
#include <iostream>
#include "scheduler.h"

namespace scheduler {


Scheduler::Scheduler(const mpi_communicator& mpi_comm, const AbstractTask& theTask)
{
  mpi_comm.recv(mpi_comm.master(),MP_task_parms,input.task_params());
  theWorker = theTask.make_worker(input.task_params());
  //theWorker->run(input.task_params());
  std::cout << "I am slave " << mpi_comm.rank() << "\n";
  //std::cout << "id = " << input.task_params().task_id() << "\n";
  //std::cout << "size = " << input.task_params().task_size() << "\n";
}

int Scheduler::run(void) 
{
  // Task task;
  // valid = input.read_params(0);
  // task.init_task_param()
  return 0;
}

} // end namespace scheduler
