/*---------------------------------------------------------------------------
* master_scheduler: Scheduler running on master
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 19:51:04
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-18 00:37:45
*----------------------------------------------------------------------------*/
#include <ctime>
#include <string>
#include <sstream>
#include "scheduler.h"

namespace scheduler {

MasterScheduler::MasterScheduler(int argc, const char *argv[], const AbstractTask& theTask)
  : Scheduler(), cmdarg(argc, argv), input(), task_size{0}
{
  if (cmdarg.valid()) {
    if (input.read_inputs(cmdarg.filename())) {
      task_size = input.task_size();
      input.init_task_parameters(parms);
      // set values for the first task to construct the worker
      input.set_task_parameters(parms, 0);
      std::cout << " starting..." << std::endl;
      theWorker = theTask.make_worker(parms);
      valid_ = true;
    }
  }
}

int MasterScheduler::run() 
{
  if (!valid_) return -1;
  auto start_time = std::chrono::steady_clock::now();
  //std::cout << " starting..." << std::endl;
  for (unsigned task=0; task<task_size; ++task) {
    //------------ timing -------------
    auto tp = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    std::string ts = std::ctime(&t);
    std::string tm = ts.substr(11,8)+" "+ts.substr(8,2)+"-"+ts.substr(4,3)+"-"+ts.substr(20,4);
    if (!cmdarg.have_option(quiet)) 
      std::cout <<" starting task "<<task+1<<" of "<< task_size<<" at "<<tm<<"\n";
    auto tstart_t = std::chrono::steady_clock::now();
    //------------------------------

    input.set_task_parameters(parms, task);
    theWorker->start(parms);
    //params << pstore(task_id);

    //-------------------------------------
    auto tend_t = std::chrono::steady_clock::now();
    if (!cmdarg.have_option(quiet)) 
      std::cout << " - finished this task in " << elapsed_time(tstart_t, tend_t) << "\n";
    //-------------------------------------
  }
  // finish
  theWorker->finish();

  auto end_time = std::chrono::steady_clock::now();
  std::cout << " finished all tasks\n"; 
  std::cout << " time: " << elapsed_time(start_time, end_time) << "\n";
  std::cout << " done!\n";
  return 0;
}

std::string MasterScheduler::elapsed_time(const std::chrono::steady_clock::time_point& start_time,
  const std::chrono::steady_clock::time_point& end_time) const
{
  std::ostringstream os;
  auto diff = end_time - start_time;
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(diff).count();
  auto days = secs/86400;
  secs = secs%86400;
  auto hrs = secs/3600;
  secs = secs%3600;
  auto mins = secs/60;
  secs = secs%60;
  if (days > 0) {
    os << days << "d:" << hrs << "h:" << mins << "m";
    return os.str();
  }
  else if (hrs > 0) {
    os << hrs << "h:" << mins << "m:" << secs << "s";
    return os.str();
  }
  else {
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    mins = milli/60000;
    milli = milli%60000;
    secs = milli/1000;
    milli = milli%1000;
    os << mins << "m:" << secs << "s:" << milli << "ms";
    return os.str();
  }
}

} // end namespace scheduler
