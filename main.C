//Author: Nestor Solalinde
//Version: 1.30


// C++ include files that we need
#include "iostream"
#include <algorithm>
#include <math.h>
#include "getpot.h"
#include "multiphase.h"
#include "perf_log.h"

int main (int argc, char** argv)
{

  GetPot input_file("./inputs/validate3.in");
  Multiphase multiphase(argc,argv);
  PerfLog perf_log;
  perf_log.push("Full Program");
  
  multiphase.initialize(input_file);
  multiphase.project_solutions();
  multiphase.initial_refinement();
  multiphase.project_solutions();
  multiphase.calculate_curvature_heaviside();
  multiphase.save_global_results();
  while(multiphase.simulating())
  {
      multiphase.increase_timestep();
      perf_log.push("Reinitialization");
      multiphase.level_reinitialization();
      perf_log.pop("Reinitialization");
      
      perf_log.push("Curvature and Heaviside");
      multiphase.calculate_curvature_heaviside();
      perf_log.pop("Curvature and Heaviside");

      perf_log.push("Navier Stokes");
      multiphase.store_ns_old_solutions();
      multiphase.solve_navier_stokes();
      perf_log.pop("Navier Stokes");
      multiphase.print_cfl();
      //multiphase.save_results_to_file();

      perf_log.push("Level Set Advection");
      multiphase.level_advection();
      perf_log.pop("Level Set Advection");
      multiphase.save_global_results();
  }
  perf_log.pop("Full Program");

  return 0;
}


