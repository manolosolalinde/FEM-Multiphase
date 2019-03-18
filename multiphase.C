//Autor: Nestor Solalinde
//Version: 1.30



#include <statistics.h>
#include <bits/stl_vector.h>

#include "multiphase.h"
//



 
Multiphase::Multiphase(int argc, char** argv): //, Mesh& meshref
  init(argc, argv),
  mesh (2),
  equation_systems (mesh),
  system_lvlset (equation_systems.add_system<TransientAdvectionDiffusionSystem> ("LevelSet")),
  system_reinit (equation_systems.add_system<TransientAdvectionDiffusionSystem> ("Reinit_LevelSet")),
  system_test (equation_systems.add_system<LinearImplicitSystem> ("TEST SYSTEM")),
  system_normalx (equation_systems.add_system<TransientAdvectionDiffusionSystem> ("normalx SYSTEM")),
  system_normaly (equation_systems.add_system<TransientAdvectionDiffusionSystem> ("normaly SYSTEM")),
  system_normalz (equation_systems.add_system<TransientAdvectionDiffusionSystem> ("normalz SYSTEM")),
  //system_stokes (equation_systems.add_system<TransientNonlinearImplicitSystem> ("Navier-Stokes")),
  system_stokes (equation_systems.add_system<TransientLinearImplicitSystem> ("Navier-Stokes")),
  system_normal (equation_systems.add_system<LinearImplicitSystem>("LevelSet Normal")),
  system_refinement (equation_systems.add_system<TransientLinearImplicitSystem>("Refinement")),
  system_corrector (equation_systems.add_system<LinearImplicitSystem>("Reinit Corrector")),
  system_heaviside(equation_systems.add_system<LinearImplicitSystem>("Heaviside Function")),
  system_kurvature(equation_systems.add_system<LinearImplicitSystem>("Kurvature")),
//  system_jump_pres(equation_systems.add_system<LinearImplicitSystem>("Jump Pres System")),
//  system_jump_vel(equation_systems.add_system<LinearImplicitSystem>("Jump Vel System")),
//  system_force(equation_systems.add_system<LinearImplicitSystem>("Force System")),
//  system_wforce(equation_systems.add_system<LinearImplicitSystem>("Whole Force System")),
  obj_interface (equation_systems,"LevelSet"),
  mesh_refinement (mesh),
  myinputfile()
  
  
{

    target_name = argv[0];

    //initiallize NULL pointers
    stokes_dp_boundary_fptr = 0;
    stokes_is_dp_boundary_fptr = 0;
    stokes_dv_boundary_fptr = 0;
    stokes_is_dv_boundary_fptr = 0;
    initial_lvlset_fptr = 0;
    exodusobj_tau = 0;
    exodusobj = 0;
    

    //Default dimension
    dim = 2;

    

    system_stokes.add_variable ("u", SECOND);//SECOND
    system_stokes.add_variable ("v", SECOND);//SECOND
    system_stokes.add_variable ("p", FIRST);
    system_stokes.attach_assemble_function (assemble_stokes);
    //system_stokes.nonlinear_solver->residual = compute_residual;
    //system_stokes.nonlinear_solver->jacobian = compute_jacobian;





    system_normal.add_variable("xnormal",SECOND);//SECOND
    system_normal.add_variable("wynormal",SECOND);//SECOND
    system_normal.add_variable("kurv",SECOND);//SECOND
    system_normal.add_variable("auxf",SECOND);//SECOND
    system_normal.attach_assemble_function(assemble_normal);


    //must be in the same space as pressure
    /*system_jump_pres.add_variable("jump_pres",FIRST);
    system_jump_pres.attach_assemble_function(assemble_jump);

    system_jump_vel.add_variable("jump_vel",SECOND);
    system_jump_vel.attach_assemble_function(assemble_jump);


    system_force.add_variable("forcex",SECOND);
    system_force.add_variable("forcey",SECOND);
    if (dim==3)
        system_force.add_variable("forcez",SECOND);
    system_force.attach_assemble_function(assemble_force);
    
    system_wforce.add_variable("wforcex",SECOND);
    system_wforce.add_variable("wforcey",SECOND);
    if (dim==3)
        system_wforce.add_variable("wforcez",SECOND);
    system_wforce.attach_assemble_function(assemble_whole_force);
    */

    


    //TODO eliminate system refinement, use grad_heaviside
    system_refinement.add_variable("refinement",SECOND);//SECOND
    system_refinement.attach_assemble_function(assemble_refinement);

    system_corrector.add_variable("lambda sussman",SECOND);//SECOND
    system_corrector.attach_assemble_function(assemble_corrector);

    system_heaviside.add_variable("H",FIRST);
    system_heaviside.attach_assemble_function(assemble_heaviside);

    system_kurvature.add_variable("k",SECOND);//SECOND
    system_kurvature.attach_assemble_function(assemble_kurvature);

    system_lvlset.add_variable("Level",SECOND);//SECOND
    system_lvlset.attach_variable_for_assemble(system_stokes,"u");
    system_lvlset.attach_variable_for_assemble(system_stokes,"v");
    system_lvlset.attach_velocity_function(lvlset_velocity);
    system_lvlset.data.theta = 0.5;

    system_reinit.add_variable("SgnLevel",SECOND);//SECOND
    system_reinit.attach_variable_for_assemble(system_lvlset,"Level");
    //system_reinit.attach_variable_for_assemble(system_reinit,"SgnLevel","old_grad_mean");
    system_reinit.attach_velocity_function(reinit_velocity);
    system_reinit.attach_force_function(reinit_force);
    system_reinit.attach_epsilon_function(reinit_epsilon);
    system_reinit.data.theta = 1;


    system_normalx.add_variable("KROBEXss",SECOND);
    system_normalx.attach_variable_for_assemble(system_reinit,"SgnLevel","grad");
    system_normalx.attach_variable_for_assemble(system_reinit,"SgnLevel");
    system_normalx.attach_velocity_function(normalx_velocity);
    system_normalx.attach_force_function(normalx_force);
    system_normalx.data.theta = 1;
    system_normalx.data.dt = 0.1;
    system_normalx.data.epsilon = 0.0001;

    system_normaly.add_variable("KROBEYss",SECOND);
    system_normaly.attach_variable_for_assemble(system_reinit,"SgnLevel","grad");
    system_normaly.attach_variable_for_assemble(system_reinit,"SgnLevel");
    system_normaly.attach_velocity_function(normaly_velocity);
    system_normaly.attach_force_function(normaly_force);
    system_normaly.data.theta = 1;
    system_normaly.data.dt = 0.1;
    system_normaly.data.epsilon = 0.0001;

    //TODO3D FIXME
    system_normalz.add_variable("KROBEZss",SECOND);


    system_test.add_variable("OMNI",SECOND);
    system_test.attach_assemble_function(assemble_test);


}

void Multiphase::initialize(GetPot& input_file){


    myinputfile.absorb(input_file);

    //MESH SETUP AND INITIALIZATION
    {
    dim = input_file("mesh_dimension",2);
    int mesh_elements_x = input_file("mesh_elements_x",9);
    int mesh_elements_y = input_file("mesh_elements_y",12);
    Real mesh_xlim_a = input_file("mesh_xlim_a",-1.5);
    Real mesh_xlim_b = input_file("mesh_xlim_b",1.5);
    Real mesh_ylim_a = input_file("mesh_ylim_a",0);
    Real mesh_ylim_b = input_file("mesh_ylim_b",4);
    MeshTools::Generation::build_square (mesh,
                                       mesh_elements_x, mesh_elements_y,
                                       mesh_xlim_a, mesh_xlim_b,
                                       mesh_ylim_a, mesh_ylim_b,
                                       QUAD9);
    equation_systems.init();
    equation_systems.print_info();
    }
//    ParmetisPartitioner partitioner;
//    partitioner.partition(mesh);
    mesh.print_info(std::cout);
    std::cout << "Number of partitions = " << mesh.n_partitions() << std::endl;


    //RESULT STORAGE PARAMETERS
    {
        save_interval = input_file("save_interval",1);
        save_tau_interval = input_file("save_tau_interval",1);
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        std::string inputfile = input_file[0];
        inputfile.erase(inputfile.begin(),inputfile.begin()+9);
        output_foldername.clear();
        output_foldername += "/home/manolo/ManoloYacyreta/outputs/output_";
        output_foldername += inputfile;
        output_foldername.erase(output_foldername.end()-3,output_foldername.end());
        output_foldername += "_";
        output_foldername += asctime(timeinfo);
        output_foldername.erase(output_foldername.end()-1,output_foldername.end());
        size_t found;
        found=output_foldername.find_first_of(" ");
        while (found!=std::string::npos)
        {
            output_foldername[found]='_';
            found=output_foldername.find_first_of(" ",found+1);
        }
        std::string auxvar = input_file("custom_output_foldername",output_foldername);
        output_foldername = auxvar;
        //output_foldername.erase(output_foldername.end()-1,output_foldername.end());
        std::string auxstr("mkdir ");
        auxstr += output_foldername;
        if (system(auxstr.c_str()) == -1) { // Create the directory
        std::cerr << "Error: " << strerror(1) << " - Cannot create directory";
        libmesh_error();
        }
        auxstr.clear();
        auxstr = "cp ";
        auxstr += input_file[0];
        auxstr += " ";
        auxstr += output_foldername;
        if (system(auxstr.c_str()) == -1) { // Create the directory
        std::cerr << "Error: " << strerror(2) << " - Cannot copy inputfile to outputs";
        libmesh_error();
        }


        std::string auxstr2 = input_file[0];
        size_t found2;
        found2 = auxstr2.find_last_of("/");
        auxstr2 = auxstr2.substr(found2+1);
        auxstr2.erase(auxstr2.end()-3,auxstr2.end());

        auxstr2.erase(auxstr2.begin(),auxstr2.end()-6);
        out_anim_filename.clear();
        out_anim_filename += output_foldername;
        out_anim_filename += "/out_anim_";
        out_anim_filename += auxstr2;
        out_anim_filename += ".ex2";
        out_refined_filename.clear();
        out_refined_filename += output_foldername;
        out_refined_filename += "/out_";
        out_refined_filename += auxstr2;
        //out_refined_filename += ".ex2";

        out_tau_filename.clear();
        out_tau_filename += output_foldername;
        out_tau_filename += "/out_anim_tau.ex2";
        exodusobj = new ExodusII_IO(mesh);
        exodusobj->write_timestep(out_anim_filename.c_str(),equation_systems,1,0);
    }


    //PHISICAL PROPERTIES
    equation_systems.parameters.set<Real>("RhoA") = input_file("RhoA",0.01);//internal relative fluid
    equation_systems.parameters.set<Real>("RhoB") = input_file("RhoB",1.00);//external relative fluid
    equation_systems.parameters.set<Real>("MuA") = input_file("MuA",0.5);//internal relative fluid
    equation_systems.parameters.set<Real>("MuB") = input_file("MuB",1.00);//external relative fluid
    Re = input_file("Re",1.0);//5680.5;//3.1944; //Adimensional Reynolds number
    We = input_file("We",1.0);//102040.8;//1.0204;//Adimensional Weber number
    Fr = input_file("Fr",1.0);//10.1;//0.319;//Adimensional Froude number
    equation_systems.parameters.set<Real>("Froude") = Fr;
    equation_systems.parameters.set<Real>("Reynolds") = Re;
    equation_systems.parameters.set<Real>("Weber") = We;
    equation_systems.parameters.set<Number>("Gravity On") = input_file("gravity_on",1);



    //CONTROL PARAMETERS
    convergence_tolerance = input_file("convergence_tolerance",4.e-9);
    nl_convergence = true;
    sim_time = 0;
    t_step = 0;
    system_lvlset.data.epsilon = input_file("artificial_diffusion_lvl",0.001);//std = 0.02
    system_reinit.data.epsilon = input_file("artificial_diffusion_reinit",0.00001);//std = 0.015
    //TODO 3D
    n_timesteps = input_file("n_timesteps",3000);
    cfl_ceil = input_file("cfl_ceil",0.8);
    cfl_bottom = input_file("cfl_bottom",0.2);
    n_nonlinear_steps_ceil = input_file("n_nonlinear_steps_ceil",6);
    n_nonlinear_steps_bottom = input_file("n_nonlinear_steps_bottom",2);
    dt = input_file("dt",0.0005); //std = 0.005
    dt_min = input_file("dt_min",0.000005);
    variable_dt = input_file("variable_dt",0); //set this =1 to use variable time steps
    equation_systems.parameters.set<Real>("CFL")=0;
    equation_systems.parameters.set<Real> ("time") = sim_time;
    equation_systems.parameters.set<Real>("stokes_tau_factor") = input_file("stokes_tau_factor",2.0);
    equation_systems.parameters.set<Real>("Artificial Diffusion reinit") = system_reinit.data.epsilon;


    //REFINEMENT
    refinement_x = input_file("refinement_x",1.2);
    int uniform_refinements = input_file("uniform_refinements",0);
    local_tolerance = input_file("local_tolerance",0.05);
    //global_refinement_tolerance = input_file("global_refinement_tolerance",1.0);
    lvlset_error_ceil = input_file("lvlset_error_ceil",2.0);
    max_r_steps = input_file("max_r_steps",2); //std=1 refinement steps per time step.
    base_max_h_level = input_file("base_max_h_level",1);//std=1
    equation_systems.parameters.set<Real>("local tolerance")=local_tolerance;
    mesh_refinement.coarsen_threshold() = input_file("coarsen_threshold",0.1);
    mesh_refinement.max_h_level() = base_max_h_level;
    equation_systems.parameters.set<unsigned int>("max_h_level") = base_max_h_level;



    //INTERFACE AND BOUNDARY DEFINITION
    {
        Real epsilon = input_file("epsilon",0.2);
        Real epsilon_min = input_file("epsilon_min",0.1);
        Real epsilon_ref = epsilon*input_file("epsilon_ref",1.25);
        obj_interface.setEpsilon(epsilon);
        obj_interface.setEpsilon_min(epsilon_min);
        obj_interface.setEpsilon_ref(epsilon_ref);
        obj_interface.set_refinement_level(base_max_h_level);
        initial_lvlset_function = input_file("initial_lvlset_function",0);
        obj_interface.print_info(std::cout);
        Real advection_velocity = input_file("test_velocity",300.00);
        equation_systems.parameters.set<bool>("Set Constant Velocity") = input_file("test_advection",0); // Standard = 0;
        equation_systems.parameters.set<RealVectorValue>("velocity") = RealVectorValue (0, advection_velocity);
        equation_systems.parameters.set<bool>("boundary_close_box") = input_file("boundary_close_box",1);
        equation_systems.parameters.set<bool>("lid_moving") = input_file("lid_moving",0);//std =0
        equation_systems.parameters.set<Real>("lid_velocity") = input_file("lid_velocity",1.0);
        equation_systems.parameters.set<Interface *>("interface_obj_ptr") = &(obj_interface);
        //stokes_dv_boundary_fptr = dirichlet_boundary1;

        unsigned int boundary_case = input_file("boundary_case",0);
        switch (boundary_case)
        {
            case 0:
                break;
            case 1:
                attach_stokes_v_dirichlet_functions(is_dirichlet_boundary_v_stokes1,dirichlet_boundary_v_stokes1);
                attach_stokes_p_dirichlet_functions(is_dirichlet_boundary_p_stokes1,dirichlet_boundary_p_stokes1);
                break;
            case 2:
                attach_stokes_v_dirichlet_functions(is_dirichlet_boundary_v_stokes2,dirichlet_boundary_v_stokes2);
                break;
            default:
                break;
        }



        equation_systems.parameters.set<RealVectorValue (*)(const Real , const Real , const Real)>("stokes_dv_boundary_fptr") = stokes_dv_boundary_fptr;
        equation_systems.parameters.set<Real (*)(const Real , const Real , const Real)>("stokes_dp_boundary_fptr") = stokes_dp_boundary_fptr;
        equation_systems.parameters.set<bool (*)(const Real , const Real , const Real)>("stokes_is_dv_boundary_fptr") = stokes_is_dv_boundary_fptr;
        equation_systems.parameters.set<bool (*)(const Real , const Real , const Real)>("stokes_is_dp_boundary_fptr") = stokes_is_dp_boundary_fptr;
        //Initialize Interface object
    }



    //REINITIALIZATION PARAMETERS
    reinit_tolerance = input_file("reinit_tolerance",2.e-3);
    initial_dt_tau = input_file("initial_dt_tau",0.1); //std = 0.1
    system_normalx.data.dt = input_file("normalx_data_dt",0.1);
    system_normaly.data.dt = input_file("normalx_data_dt",0.1);;
    //TODO 3D
    n_initial_tau_steps = input_file("n_initial_tau_steps",10);//std =80
    n_tau_steps = input_file("n_tau_steps",1); //std = 1
    reinit_record_step = input_file("reinit_record_step",n_timesteps);//std= n_timesteps
    last_n_tau_steps = input_file("last_n_tau_steps",60);
    sussman_corrector = input_file("sussman_corrector",1);
    bool ksupg_reinit_on = input_file("ksupg_reinit_on",1);
    equation_systems.parameters.set<Real>("velocity_low_limit") = input_file("velocity_low_limit",0.0);
    equation_systems.parameters.set<Real>("diffusion_factor_lowveloc") = input_file("diffusion_factor_lowveloc",4.0);
    equation_systems.parameters.set<Real> ("dt_tau") = initial_dt_tau;
    equation_systems.parameters.set<Real> ("dt")   = dt;
    equation_systems.parameters.set<bool> ("ksupg_reinit_on") = ksupg_reinit_on;




    //NONLINEAR SOLVER
    n_nonlinear_steps = input_file("n_nonlinear_steps",50);
    l_step = n_nonlinear_steps_ceil; //nonlinear steps counter
    nonlinear_absolute_tolerance = input_file("nonlinear_absolute_tolerance",1.e-8);
    nonlinear_relative_tolerance = input_file("nonlinear_relative_tolerance",1.e-12);
    nonlinear_step_tolerance = input_file("nonlinear_step_tolerance",1.e-6);
    linsolver_max_iterations = input_file("linsolver_max_iterations",2000);
    equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = linsolver_max_iterations;
    equation_systems.parameters.set<unsigned int> ("nonlinear solver maximum iterations") = n_nonlinear_steps;
    equation_systems.parameters.set<Real>("maxtsupg")=1e10;
//    #ifdef LIBMESH_HAVE_PETSC
//        std::stringstream atol;
//        atol << nonlinear_absolute_tolerance;//||F|| < atol
//        std::string atol2 = atol.str();
//        std::stringstream stol;            // PNORM_RELATIVE
//        stol << nonlinear_step_tolerance; // Newton computed step size small; || delta x || < stol
//        std::string stol2 = stol.str();
//        std::stringstream rtol;
//        rtol << nonlinear_relative_tolerance;//||F|| < rtol*||F_initial||
//        std::string rtol2 = rtol.str();
//        PetscOptionsSetValue("-snes_stol",stol2.c_str());
//        PetscOptionsSetValue("-snes_atol",atol2.c_str());
//        PetscOptionsSetValue("-snes_rtol",rtol2.c_str());
//    #endif



    //EQ.SYSTEM REINITIALIZATION AND UNIFORM REFINEMENT
    mesh_refinement.uniformly_refine(uniform_refinements);
    equation_systems.reinit();
    

}

void Multiphase::project_solutions(){

    //no borrar este comentario, se queda como cicatriz
    //    Number (Multiphase::*fptr)(const Point& p,const Parameters& parameters,
    //            const std::string& sys_name,
    //            const std::string& unknown_name);
    //    fptr = &Multiphase::exact_value;

    if (initial_lvlset_fptr == NULL)
        initial_lvlset_fptr = function_one_bubble;

    Number (*fptr)(const Point&, const Parameters&,
            const std::string&,const std::string&);
    switch (initial_lvlset_function)
    {
        case 0:
            fptr = initial_lvlset_fptr;
            break;
        case 1:
            fptr = function_two_bubbles;
            break;
        case 2:
            fptr = function_testreinit;
            break;
        case 3:
            fptr = function_eliptic_bubble;
            break;
        case 4:
            fptr = function_one_bubble;
            break;
        case 5:
            fptr = function_one_phase;
            break;
        default:
            fptr = initial_lvlset_fptr;
    }

    //Inicializacion de variables
    system_reinit.project_solution(fptr, NULL, equation_systems.parameters);
    system_lvlset.project_solution(fptr, NULL, equation_systems.parameters);
    
    system_reinit.reinit();
    system_lvlset.reinit();
    equation_systems.reinit();

}

void Multiphase::initial_refinement(){

    lvlset_error.clear();
    ErrorVector heaviside_error;
    heaviside_error.clear();
    
    KellyErrorEstimator error_estimator;
    error_estimator.error_norm = L2;
    Real initial_ref_rel_tolerance = myinputfile("initial_ref_rel_tolerance",0.5);
    local_tolerance = local_tolerance*initial_ref_rel_tolerance;
    bool keep_refining = true;
    for (int i=0;i<15;i++)
    {
	 
        //level_reinitialization(1);
        system_heaviside.solve();
	 
	 
        error_estimator.estimate_error (system_heaviside, lvlset_error);
 

        Real sqrt_n = std::sqrt(static_cast<Real>(mesh.n_active_elem()));
        

          std::cout << "local_tolerance = " << local_tolerance
            << "   lvlset_error*refinement_x = " << lvlset_error.l2_norm()/sqrt_n*refinement_x;

          if (local_tolerance < lvlset_error.l2_norm()/sqrt_n*refinement_x)
          {
              std::cout << "   Error is not acceptable\n";
              keep_refining = refinement(1,false,true);
              if (keep_refining==false)
              {
                  std::cout << "   Error is acceptable\n";
                  break;
              }

          }
          else
          {
              std::cout << "   Error is acceptable\n";
              break;
          }

    }
    local_tolerance = local_tolerance/initial_ref_rel_tolerance;
}

bool Multiphase::refinement(int ref_steps,bool refine_by_stokes_error,bool refine_interface){ ///recien copiado, arreglar


    Real worst_error = 0;
    bool keep_refining = true;
    for (unsigned int r_step=0; r_step<ref_steps && keep_refining; r_step++)
    {
      std::cout << "  Refining the mesh...";
      Real global_tolerance = local_tolerance*std::sqrt(static_cast<Real>(mesh.n_active_elem()));
      mesh_refinement.absolute_global_tolerance() = global_tolerance;
    
      {


          if(refine_interface)
          {
              //TODO
          }
        
          if(refine_by_stokes_error)
          {
              if ((stokes_error.size() != mesh.n_active_elem() ) || r_step>0)
              {
              stokes_error.clear();
              stokes_error = calculate_stokes_error();
              }
          }

          if(!refine_by_stokes_error && !refine_interface)
          {
              libmesh_error();
          }
          

          ErrorVector total_error;
          if (refine_by_stokes_error && refine_interface)
          {
              total_error.resize(std::max(lvlset_error.size(),stokes_error.size()));
              for (unsigned int i=0; i<lvlset_error.size();i++)
              {
                  total_error[i]=std::max(lvlset_error[i],stokes_error[i]);
              }
          }
          else
          {
              if(refine_by_stokes_error)
                  total_error = stokes_error;
              if(refine_interface)
                  total_error = lvlset_error;
          }
          
          worst_error = total_error.maximum();
          std::cout << " Worst element error = " << worst_error
                          << ", mean = " << total_error.mean();

          mesh_refinement.flag_elements_by_error_tolerance(total_error);
          keep_refining = mesh_refinement.refine_and_coarsen_elements();
          equation_systems.reinit ();

                     std::cout << "...mesh now have " << mesh.n_active_elem()
                        << " active elements" <<std::endl;
          
       }
    }
    if (worst_error<=local_tolerance || keep_refining==false)
    {
        keep_refining=false;
    }
    else
    {
        keep_refining=true;
    }

    return keep_refining;

}

bool Multiphase::refinement(void){
    bool keep_refining = refinement(max_r_steps);
    return keep_refining;

}

void Multiphase::level_reinitialization(int numberOfSteps){ 

            AutoPtr<NumericVector<Number> > last_reinit_soln (system_reinit.solution->clone());
            AutoPtr<NumericVector<Number> > last_corrector_soln (system_corrector.solution->clone());

            // pass solution of lvlset advection to reinitialization funtion
            *system_reinit.current_local_solution = *system_lvlset.current_local_solution;

            Real dt_tau = initial_dt_tau;
            Real tau_time=0;
            Real norm_reinit_delta = 1.e3;
            Real old_norm_reinit = 1.e4;

            if ((t_step)==reinit_record_step)
            {
                exodusobj_tau = new ExodusII_IO(mesh);
                numberOfSteps = last_n_tau_steps;
                exodusobj_tau->write_timestep(out_tau_filename.c_str(),equation_systems,1,tau_time);
            }
            unsigned int tau_step;

        for (tau_step = 0; tau_step <numberOfSteps; tau_step++)
        {
            last_reinit_soln->zero();
            last_reinit_soln->add(*system_reinit.solution);
            last_corrector_soln->zero();
            last_corrector_soln->add(*system_corrector.solution);

            if (norm_reinit_delta >= 0.8)
            {
                dt_tau = dt_tau*1;
                //std::cout << "Reducing dt_tau = " << dt_tau << std::endl;
            }

            tau_time += dt_tau;
            equation_systems.parameters.set<Real> ("dt_tau") = dt_tau;
            equation_systems.parameters.set<Real> ("tau_time") = tau_time;
            system_reinit.data.dt = dt_tau;
            system_reinit.data.time = tau_time;

            *system_reinit.old_local_solution = *system_reinit.current_local_solution;

            system_reinit.solve();
            
            system_corrector.solve();
            Real max_corrector = (system_corrector.current_local_solution)->linfty_norm();

            //TODO why (sussman_corrector == 1) && (max_corrector<0.15*obj_interface.getDeltax_min())
            if (sussman_corrector == 1)//check if we are using sussman corrector or not
            {
                *system_reinit.solution += *system_corrector.current_local_solution;
                /*MeshBase::const_node_iterator it = mesh.active_nodes_begin();
                const MeshBase::const_node_iterator it_end = mesh.active_nodes_end(); 
                for ( ; it != it_end;++it)
                {
                    Real reinit_solution = (system_reinit.solution)->operator ()((*it)->id());
                    Real corrector_solution = (system_corrector.solution)->operator ()((*it)->id());
                    (system_reinit.solution)->set((*it)->id(),reinit_solution+corrector_solution);
                }*/

            }
            //recalculate reinit according to sussman
            //correct_reinit (equation_systems);
            //copy solution to current_local_solution
            
            
            system_reinit.update();

            if((t_step)==reinit_record_step && (tau_step+1)%save_tau_interval==0)
                exodusobj_tau->write_timestep(out_tau_filename.c_str(),equation_systems,(tau_step+1)/save_tau_interval+1,tau_time);

            //check convergence
            last_reinit_soln->add (-1., *system_reinit.solution);
            last_reinit_soln->close();
            last_corrector_soln->add (-1., *system_corrector.solution);
            last_corrector_soln->close();
            old_norm_reinit = norm_reinit_delta;
            norm_reinit_delta = (last_reinit_soln->l2_norm())/(*system_reinit.solution).l2_norm();


            std::cout << " Solving step " << tau_step+1 <<" of system_reinit... "
                      << " Reinit convergence |sgn - sgn_old|/|sgn| = "
                      << norm_reinit_delta << " dt_tau = " << dt_tau << std::endl;

            std::cout << "      |lambda|inf/dt = "
                      << max_corrector/dt_tau
                      << "      |lambda - lambdaold|/dt = "
                      << last_corrector_soln->linfty_norm()/dt_tau
                      << std::endl;
            if (norm_reinit_delta < reinit_tolerance)
            {
              std::cout << "   Reinitialization function converged at step "
                        << tau_step+1
                        << std::endl;
              break;
            }

        }

        //SETTING LVLSET OLD LOCAL SOLUTION = REINIT CURRENT LOCAL SOLUTION (THETA=1)
        // Note that reinit will change for the curvature calculation
        *system_lvlset.old_local_solution = *system_reinit.current_local_solution;

            
        if ((t_step)==reinit_record_step)
        {
            delete exodusobj_tau;
        }

}

void Multiphase::level_reinitialization(void){ 
    if (t_step==1)
        level_reinitialization(n_initial_tau_steps);
    else
        level_reinitialization(n_tau_steps);
}

void Multiphase::calculate_curvature_heaviside(){

        //Here we are reinitializing reinit without correction. The effect is
        //not cumulative since we are advecting the corrected solution

        Real myepsilon = system_reinit.data.epsilon;

        //system_reinit.data.epsilon = 0.005;
        //system_reinit.data.dt = 2.0;
        system_normalx.data.epsilon = system_reinit.data.epsilon;
        system_normaly.data.epsilon = system_reinit.data.epsilon;
        //system_normalx.data.dt = 0.01;//specified in initialization
        //system_normaly.data.dt = 0.01;

        //FIXME este reinit no puede ser activado, repercute sobre la funcion de
        // reinicializacion original y no se porque
        //system_reinit.solve();
        
        system_normalx.solve();
        system_normaly.solve();
        system_kurvature.solve();

        system_reinit.data.epsilon = myepsilon;
        
        system_normal.solve();//FIXME - WE ARE NOT USING THIS SYSTEM
        
        system_heaviside.solve();

        //system_wforce.solve();

        //system_jump_pres.solve();
        //system_jump_vel.solve();
        //system_force.solve();
        

}

bool Multiphase::simulating(){

    
    Real norm_rel=5000;
    if (t_step>0)
    {
    AutoPtr<NumericVector<Number> > last_nonlinear_soln (system_stokes.current_local_solution->clone());
          last_nonlinear_soln->zero();
          last_nonlinear_soln->add(*system_stokes.current_local_solution);
          last_nonlinear_soln->add (-1., *system_stokes.old_local_solution);
          Real norm_diff = last_nonlinear_soln->l2_norm();
          norm_rel = norm_diff/(system_stokes.current_local_solution)->l2_norm();
          std::cout << " Now ||u - u_old||/||u|| = "<<norm_rel<<std::endl;
          std::cout << "     ||u - u_old|| = "<< norm_diff <<std::endl;
	   if (l_step == n_nonlinear_steps)
		norm_rel = 1e10;

        if ((t_step < n_timesteps && norm_rel > convergence_tolerance*dt) || norm_diff==0)
        {
            return true;
        }
        else
        {
            return false;
        }

    }
    else
    {
        return true; //case t_step=0

    }

    return true;
    
}

void Multiphase::set_n_timesteps(unsigned int numberOfTimesteps){
    n_timesteps = numberOfTimesteps;
    reinit_record_step = numberOfTimesteps;
}

int Multiphase::get_n_timesteps(){
    return n_timesteps;
}

void Multiphase::stop(){
    t_step = n_timesteps;//last timestep
}

void Multiphase::increase_timestep(){
    //Incrementing timestep counter
    t_step++;
    calculate_timestep();
}

void Multiphase::calculate_timestep(){


    myinputfile.parse_input_file(myinputfile[0]);
        if(myinputfile.size()!=0)
        {
            cfl_ceil = myinputfile("cfl_ceil",0.8);
            cfl_bottom = myinputfile("cfl_bottom",0.2);
            n_nonlinear_steps_ceil = myinputfile("n_nonlinear_steps_ceil",6);
            n_nonlinear_steps_bottom = myinputfile("n_nonlinear_steps_bottom",2);
        }


            Real cfl = equation_systems.parameters.get<Real>("CFL");
            //decreasing time step if needed
            if ((cfl > cfl_ceil || l_step>n_nonlinear_steps_ceil) && variable_dt==1 && dt>=dt_min)//l_step+1 >= l_step_ceil
            {

                dt = dt*0.80;
                std::cout << " Reducing time interval dt to: " << dt << std::endl;
            }
            //incrementing time step if needed
            else
            if ( (cfl < cfl_bottom && l_step<n_nonlinear_steps_bottom) && variable_dt==1)//l_step+1 <= l_step_bottom
            {
                dt = dt/0.80;
                std::cout << " Increasing time interval dt to: " << dt << std::endl;
            }

            sim_time += dt; 
            equation_systems.parameters.set<Real> ("time") = sim_time;
            equation_systems.parameters.set<Real> ("dt")   = dt;
            std::cout << "\nSolving time step ";
            {
                OStringStream out;
                OSSInt(out,2,t_step);
                out << ", time=";
                OSSRealzeroleft(out,6,3,sim_time);
                out <<  "...";
                std::cout << out.str() << "with dt: " << dt << std::endl;
            }

}

/*void Multiphase::solve_navier_stokes(){

      const Real initial_linear_solver_tol = 1.e-6;
      equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;

      nl_convergence = false;
      

      for (unsigned int l = 0; l<max_r_steps;l++)
      {
          std::cout << " Solving Navier-Stokes nonlinear system...\n";
          system_stokes.solve();
          unsigned int auxl_step;
          auxl_step = system_stokes.n_nonlinear_iterations();

          if(l==0)
              l_step=auxl_step;

          nl_convergence = system_stokes.nonlinear_solver->converged;

          if (nl_convergence == true)
                std::cout << "Nonlinear solver Converged after " << auxl_step << " iterations\n";
          else
          {
              std::cout << "Nonlinear solver did not converged\n";
              l_step = n_nonlinear_steps;
          }

          Real sqrt_n = std::sqrt(static_cast<Real>(mesh.n_active_elem()));
          stokes_error = calculate_stokes_error();
          //stokes_error.plot_histogram("histogram.txt",20);


          std::cout << "local_tolerance = " << local_tolerance
                    << "   Error_mean*refinement_x = " << stokes_error.l2_norm()/sqrt_n*refinement_x;
          
          if (local_tolerance < stokes_error.l2_norm()/sqrt_n*refinement_x || nl_convergence==false)
          {
              std::cout << "   Error is not acceptable\n";
              refinement(1);
          }
          else
          {
              std::cout << "   Error is acceptable\n";
              break;
          }
      }

      std::cout << "  Max tsupg = " << equation_systems.parameters.get<Real>("maxtsupg") << std::endl;


}
*/

void Multiphase::solve_navier_stokes(){

      const Real initial_linear_solver_tol = 1.e-6;
      equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;

      bool nl_convergence = false;

          //calculate_timestep();
          //nl_convergence = multiphase.solve_navier_stokes();

          // Update the nonlinear solution.

      int a_step=0;

      // Now we begin the nonlinear loop
      for (l_step=0; l_step<n_nonlinear_steps; ++l_step)
        {
         AutoPtr<NumericVector<Number> > last_nonlinear_soln (system_stokes.solution->clone());
          last_nonlinear_soln->zero();
          last_nonlinear_soln->add(*system_stokes.solution);

          // Assemble & solve the linear system.
          system_stokes.solve();

          // Compute the difference between this solution and the last
          // nonlinear iterate.
          last_nonlinear_soln->add (-1., *system_stokes.solution);

          // Close the vector before computing its norm
          last_nonlinear_soln->close();


          // Compute the l2 norm of the difference
          const Real norm_delta = last_nonlinear_soln->l2_norm();
          const Real norm_u = (*system_stokes.solution).l2_norm();


          // How many iterations were required to solve the linear system?
          const unsigned int n_linear_iterations = system_stokes.n_linear_iterations();


          // What was the final residual of the linear system?
          const Real final_linear_residual = system_stokes.final_linear_residual();

          // Print out convergence information for the linear and
          // nonlinear iterations.
          std::cout << "Linear solver converged at step: "
                    << n_linear_iterations
                    << ", final residual: "
                    << final_linear_residual
                    << "  Nonlinear convergence: ||u - u_old||/||u|| = "
                    << norm_delta/norm_u << "   ||u|| = " << norm_u
                    << std::endl;


          if ((norm_delta/norm_u < nonlinear_step_tolerance) &&
              (system_stokes.final_linear_residual() < initial_linear_solver_tol))
            {
              std::cout << " Nonlinear solver converged at step "
                        << l_step
                        << std::endl;
              nl_convergence = true;
              break;
            }

            equation_systems.parameters.set<Real> ("linear solver tolerance") =
            //std::min(final_linear_residual/100, initial_linear_solver_tol);
            std::min(Utility::pow<2>(final_linear_residual), initial_linear_solver_tol);
//            refinement(1);

    //system_test.solve();//TODO to be deleted


       }//closes nonlinear loop

          Real sqrt_n = std::sqrt(static_cast<Real>(mesh.n_active_elem()));
          stokes_error = calculate_stokes_error();
          //stokes_error.plot_histogram("histogram.txt",20);


          std::cout << "local_tolerance = " << local_tolerance
                    << "   Error_mean*refinement_x = " << stokes_error.l2_norm()/sqrt_n*refinement_x;

          if (local_tolerance < stokes_error.l2_norm()/sqrt_n*refinement_x || nl_convergence==false)
          {
              std::cout << "   Error is not acceptable\n";
              refinement(1);
          }
          else
          {
              std::cout << "   Error is acceptable\n";
              //break;
          }   

      std::cout << "Max tsupg = " << equation_systems.parameters.get<Real>("maxtsupg") << std::endl;



//      if (nl_convergence == false)
//          stop();

}

ErrorVector Multiphase::calculate_stokes_error(){


    
          ErrorVector stokes_error;

          AutoPtr<ErrorEstimator> error_estimator;
          //UniformRefinementEstimator *u = new UniformRefinementEstimator;
          KellyErrorEstimator *u = new KellyErrorEstimator;
          u->error_norm = L2;
          error_estimator.reset(u);

          Real u_weight = myinputfile("weight_u",1.0);
          Real p_weight = myinputfile("weight_p",1.0);

          // Calculate error based on u and v (and w?) but not p
          std::vector<Real> weights(2,u_weight);  // u, v
          if (2 == 3)
            weights.push_back(1.0);          // w
          weights.push_back(p_weight);            // p
          // Keep the same default norm type.
          std::vector<FEMNormType>
            norms(1, error_estimator->error_norm.type(0));
          error_estimator->error_norm = SystemNorm(norms, weights);

          error_estimator->estimate_error(system_stokes, stokes_error);


          return stokes_error;


}

void Multiphase::print_cfl(){

          //last_nonlinear_soln->zero();
          //last_nonlinear_soln->add(*system_stokes.old_local_solution);
          AutoPtr<NumericVector<Number> > last_nonlinear_soln (system_stokes.old_local_solution->clone());
          last_nonlinear_soln->add (-1., *system_stokes.current_local_solution);
          last_nonlinear_soln->close();
          Real norm_delta = last_nonlinear_soln->l2_norm();
          //std::cout<< "||stokes.old_local_solution - stokes.current_local_solution|| = "
          //        << norm_delta << std::endl;
          //norm_delta_0 = norm_delta;

          std::cout<< "  Max CFL = " <<  equation_systems.parameters.get<Real>("CFL") << std::endl;
          std::cout<< "  Max solution value (u,v,p) = " << system_stokes.current_local_solution->max() <<std::endl;
}

void Multiphase::store_ns_old_solutions(){
        //stores the old solutions for the navier stokes equation
        *system_stokes.older_local_solution = *system_stokes.old_local_solution;
        *system_stokes.old_local_solution = *system_stokes.current_local_solution;
}

void Multiphase::level_advection(){
        //Passing value of the lvlset function reinitialized
        //FIXME this is repeated
        *system_lvlset.current_local_solution = *system_reinit.current_local_solution;
        *system_lvlset.old_local_solution = *system_reinit.current_local_solution;
        std::cout << " Solving system_lvlset..." << std::endl;
        std::cout<< "lvlset max = " <<(*system_lvlset.old_local_solution).linfty_norm();
        system_lvlset.solve();
        std::cout<< "lvlset new max = " <<(*system_lvlset.current_local_solution).linfty_norm();
}

void Multiphase::save_global_results(){

    if(t_step==0 && max_r_steps==0)
    {
    //Initialize Exodus object
    exodusobj = new ExodusII_IO(mesh);
    exodusobj->write_timestep(out_anim_filename.c_str(),equation_systems,1,0);
    }

    //Update outfile.txt
    std::string auxstr("cp ./outfile_sim");
    auxstr += target_name[5];
    auxstr += ".txt ";
    auxstr += output_foldername;
    if (system(auxstr.c_str()) == -1) { // Copy the file
    std::cerr << "Error: " << strerror(1) << " - Cannot copy the file";
    libmesh_error();
    }

    if ( (t_step)%save_interval == 0)
  	  {
  	    OStringStream file_name2;
            file_name2 << out_refined_filename;
            file_name2 << ".";
  	    OSSRealzeroright(file_name2,4,0,(t_step)/save_interval+1);
  	    file_name2 << ".ex2";

            exodusobj->write_timestep(out_anim_filename,equation_systems,(t_step)/save_interval+1,sim_time);

  	    //GMVIO(mesh).write_equation_systems (file_name1.str(),
            //				equation_systems);
            ExodusII_IO(mesh).write_equation_systems (file_name2.str(), equation_systems);
            //ExodusII_IO (mesh).write_timestep("out_anim.ex2",equation_systems,t_step+1,time);

  	  }
}

void Multiphase::save_results_to_file(){
    
    Number alpha = get_value_on_coordinate(0,0.5,0,system_kurvature,"k");
    Number beta = get_value_on_coordinate(0,0.5,0,system_stokes,"v");
    std::cout << " Kurvature = "<< alpha <<std::endl;
    std::cout << " Velocity = "<< beta <<std::endl;
    
}

Number Multiphase::get_value_on_coordinate(Real x, Real y, Real z,const System& system, const std::string& unknown_name)
{
    AutoPtr<MeshFunction> values;
    values = AutoPtr<MeshFunction>
       (new MeshFunction(equation_systems,
                         *(system.solution),
                         system.get_dof_map(),
                         system.variable_number(unknown_name)));
    values->init();
    Point p0(x,y);
    return values->operator()(p0);

}



