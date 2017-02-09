#include <iostream>
#include "simulation.h"

int main(int argc, const char *argv[])
{
  try {
	 return scheduler::start(argc, argv, scheduler::Task<mc::VMC>());
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    return -1;
  }
  catch (...) {
    std::cout << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
}
