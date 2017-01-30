/*---------------------------------------------------------------------------
* task.h: Class representing a task to be run.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 10:57:03
* Last Modified by:   amedhi
* Last Modified time: 2015-09-28 22:50:08
*----------------------------------------------------------------------------*/
#include <iostream>
#include "worker.h"

#ifndef SCHEDULER_TASK_H
#define SCHEDULER_TASK_H

namespace scheduler {

class AbstractTask 
{
public:
  ~AbstractTask() {}
  virtual Worker* make_worker(input::Parameters& p) const = 0;
  virtual void print_copyright(std::ostream& out) const = 0;
};

template <class WORKER>
class Task : public AbstractTask
{
public:
  Task() {}
  ~Task() {}
  //virtual WORKER* make_worker(input::Parameters& p) const;
  Worker* make_worker(input::Parameters& p) const override
  {
    return new WORKER(p);
  }
  void print_copyright(std::ostream& out) const override
  {
    WORKER::print_copyright(out);
  }
};


} // end namespace input

#endif
