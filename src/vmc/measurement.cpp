/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-17 23:30:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-17 23:34:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include "simulator.h"

namespace vmc {

int Simulator::do_measurements(void)
{
  if (observables_.energy()) {
    //double e = energy_terms.sum();
    //observables.energy() << e;
  }
  return 0;
}


} // end namespace vmc
