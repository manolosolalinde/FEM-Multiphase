//Author: Nestor Solalinde
//Version: 1.30


// C++ include files that we need
#include "iostream"
#include <algorithm>
#include <math.h>
#include "getpot.h"
#include "multiphase.h"

int main (int argc, char** argv)
{

  GetPot input_file("./inputs/validate2.in");
  Multiphase multiphase(argc,argv);

  multiphase.initialize(input_file);
  //multiphase.set_n_timesteps(300);
  multiphase.project_solutions();
  multiphase.refinement(5);
  multiphase.project_solutions();
  //multiphase.level_reinitialization();
  multiphase.save_global_results();
  while(multiphase.simulating())
  {
      multiphase.increase_timestep();
      multiphase.refinement();
      //multiphase.level_reinitialization();
      multiphase.calculate_density_viscosity();
      multiphase.calculate_curvature_heaviside();
    
      multiphase.solve_navier_stokes();
      
      multiphase.print_cfl();
      multiphase.store_ns_old_solutions();
      multiphase.level_advection();
      
      multiphase.save_global_results();
      multiphase.save_results_to_file();//experimental, name to change
  }

  return 0;
}


