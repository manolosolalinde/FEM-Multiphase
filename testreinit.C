//Autor: Nestor Solalinde
//Version: 1.10


// C++ include files that we need
#include "iostream"
#include <algorithm>
#include <math.h>
#include "getpot.h"
#include "multiphase.h"
#include "perf_log.h"


int main (int argc, char** argv)
{

    PerfLog mainperflog;
    
  GetPot input_file("./inputs/testreinit.in");
  Multiphase multiphase(argc,argv);

  multiphase.initialize(input_file);
  multiphase.project_solutions();
//  multiphase.refinement(5,false,true);
//  multiphase.project_solutions();
  mainperflog.start_event("event1");
  multiphase.level_reinitialization();
  mainperflog.stop_event("event1");
  multiphase.save_global_results();
  


  return 0;
}


