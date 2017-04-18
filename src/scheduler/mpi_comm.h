/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-04-17 23:13:06
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-18 15:24:27
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

const int MP_task_parms = 0;

#ifdef HAVE_BOOST_MPI
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



} // end namespace input

#endif
