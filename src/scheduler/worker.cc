/*---------------------------------------------------------------------------
* task.cc: Implementation of the task classes.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 11:01:37
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:02:47
*----------------------------------------------------------------------------*/
#include "worker.h"

namespace scheduler {

Worker::Worker() : finished_(false), started_(false)
{
}

Worker::~Worker()
{
}

int Worker::start(input::Parameters& p)
{
  started_ = true;
  return 0;
}

void Worker::run()
{
  //if(started() && !finished_)
  dostep();
}

void Worker::halt()
{
  started_=false;
}

/*
void Task::finish()
{
  finished_=true;
}

// halt all active runs
*/


} // end namespace input
