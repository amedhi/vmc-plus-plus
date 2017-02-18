/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-17 23:30:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-18 07:59:55
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include "simulator.h"

namespace vmc {

int Simulator::do_measurements(void)
{
  if (observables_.energy()) {
    // bond energies
    /*using qn_id = model::qn_id;
    for (auto it=bond_terms_.cbegin(); it!= bond_terms_.cend(); ++it) {
      switch (it->qn_operator().id()) {
        case qn_id::cdagc_sigma:
          //term = op_upspin_hop(); 
          break;
      }
    }*/
    //double e = energy_terms.sum();
    //observables.energy() << e;
  }
  return 0;
}


} // end namespace vmc
