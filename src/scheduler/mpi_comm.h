/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-04-17 23:13:06
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-21 00:12:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MPI_COMM_H
#define MPI_COMM_H

#include <iostream>
#ifdef HAVE_BOOST_MPI
  #include <boost/mpi/environment.hpp>
  #include <boost/mpi/communicator.hpp>
#endif

namespace scheduler {

enum {MP_make_task, MP_task_params, MP_run_task, MP_task_finished,
  MP_quit_tasks};

//const int MP_task_params = 0;
//const int MP_run_task = 1;
//const int MP_quit_tasks = 2;

#ifdef HAVE_BOOST_MPI

using mpi_status = boost::mpi::status;
//  using mpi_environment = boost::mpi::environment;
//  using mpi_communicator = boost::mpi::communicator;
class mpi_environment : public boost::mpi::environment
{
public:
  using environment::environment;
  ~mpi_environment() {}
};

class mpi_communicator : public boost::mpi::communicator
{
public:
  using communicator::communicator;
  ~mpi_communicator() {}
  bool is_master(void) const { return rank()==0; }
  int master(void) const { return 0; }
private:
};

#else

class mpi_environment
{
public:
  mpi_environment() {}
  ~mpi_environment() {}
};

class mpi_communicator
{
public:
  mpi_communicator() {}
  ~mpi_communicator() {}
  bool is_master(void) const { return rank()==0; }
  int master(void) const { return 0; }
private:
  int rank(void) const { return 0; }
  int size(void) const { return 0; }
private:
};

#endif



} // end namespace scheduler

#endif
