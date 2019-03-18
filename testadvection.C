//Autor: Nestor Solalinde
//Version: 1.10


// C++ include files that we need
#include "iostream"
#include <algorithm>
#include <math.h>
#include "getpot.h"
#include "multiphase.h"


int main (int argc, char** argv)
{

  GetPot input_file("./inputs/testadvection.in");
  Multiphase multiphase(argc,argv);

  multiphase.initialize(input_file);
  multiphase.set_n_timesteps(50);
  multiphase.refinement();
  multiphase.level_reinitialization(5);
  while(multiphase.simulating())
  {
      multiphase.calculate_timestep();
      multiphase.level_advection();
      //multiphase.refinement();
      multiphase.level_reinitialization(0);
      multiphase.save_global_results();
      multiphase.increase_timestep();
      multiphase.save_results_to_file();
  }

  return 0;
}


