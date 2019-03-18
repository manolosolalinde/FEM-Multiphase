//Autor: Nestor Solalinde
//Version: 1.30


// C++ include files that we need
#include "iostream"
#include <algorithm>
#include <math.h> 
#include <time.h>

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"
#include "exact_solution.h"
//#include "mesh_data.h"
#include "perf_log.h"
#include "boundary_info.h"
#include "utility.h"


// Define the Finite Element object.
//#include "fe.h"

// Define Gauss quadrature rules.
//#include "quadrature_gauss.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"


// Define useful datatypes for finite element
// matrix and vector components.
//#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "mesh_refinement.h"
#include "error_vector.h"
//#include "kelly_error_estimator.h"
//#include "interface_refinement_estimator.h"
#include "exact_error_estimator.h"


#include "elem.h"
//#include "dof_map.h"
#include "transient_system.h"
#include "exodusII_io.h"
#include "interface.h"
#include "assemble_systems.h"
#include "getpot.h"
#include "mesh_function.h"

#include <sys/stat.h>

#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "uniform_refinement_estimator.h"

#include "myfunctions.h"

#ifdef LIBMESH_HAVE_PETSC
#include <petsc.h>
#endif

//#include "parmetis_partitioner.h"


class Multiphase {


  LibMeshInit init;

  GetPot myinputfile;

  unsigned int dim; 

  std::string target_name;
  std::string output_foldername;
  std::string out_anim_filename;
  std::string out_refined_filename;
  std::string out_tau_filename;

  Mesh mesh;

  ExodusII_IO * exodusobj_tau;
  ExodusII_IO * exodusobj;

  EquationSystems equation_systems;
  TransientAdvectionDiffusionSystem& system_lvlset;
  TransientAdvectionDiffusionSystem& system_reinit;
  LinearImplicitSystem& system_test; //TODO delete later
  TransientAdvectionDiffusionSystem& system_normalx;
  TransientAdvectionDiffusionSystem& system_normaly;
  TransientAdvectionDiffusionSystem& system_normalz;
  //TransientNonlinearImplicitSystem& system_stokes;
  TransientLinearImplicitSystem& system_stokes;
  LinearImplicitSystem & system_normal;
//  LinearImplicitSystem & system_jump_vel;
//  LinearImplicitSystem & system_jump_pres;
//  LinearImplicitSystem & system_force;
//  LinearImplicitSystem & system_wforce;
  TransientLinearImplicitSystem & system_refinement;
  LinearImplicitSystem & system_corrector;
  LinearImplicitSystem & system_heaviside;
  LinearImplicitSystem & system_kurvature;

  Interface obj_interface;
  MeshRefinement mesh_refinement;

  int initial_lvlset_function;
  Real sim_time;
  Real Re;//5680.5;//3.1944; //Adimensional Reynolds number
  Real We;//102040.8;//1.0204;//Adimensional Weber number
  Real Fr;//10.1;//0.319;//Adimensional Froude number
  int save_interval;
  int save_tau_interval;
  bool restriction; //Set =1 to enable restrictions on reinit function

  //timestep
  unsigned int n_timesteps;
  Real cfl_ceil;
  Real cfl_bottom;
  Real dt; //std = 0.005
  Real dt_min;
  bool variable_dt; //set this =1 to use variable time steps
  unsigned int t_step;
  Real convergence_tolerance;


  //Refinement
  unsigned int max_r_steps; //std=1 refinement steps per time step.
  unsigned int base_max_h_level;//std=1
  Real local_tolerance;
  Real lvlset_error_ceil;
  Real refinement_x;
  //Real global_refinement_tolerance;

  //Reinitialization
  Real initial_dt_tau; //std = 0.1
  Real reinit_record_step;//std= n_timesteps
  int n_initial_tau_steps;//std =80
  int n_tau_steps; //std = 1
  int last_n_tau_steps;
  bool sussman_corrector;

  //nonlinear solver
  bool nl_convergence;
  unsigned int n_nonlinear_steps;
  unsigned int n_nonlinear_steps_ceil;
  unsigned int n_nonlinear_steps_bottom;
  Real reinit_tolerance;
  unsigned int linsolver_max_iterations;
  unsigned int l_step;//nonlinear steps counter
  Real nonlinear_absolute_tolerance;
  Real nonlinear_relative_tolerance;
  Real nonlinear_step_tolerance;

  

  ErrorVector stokes_error;
  ErrorVector lvlset_error;

  RealVectorValue (*stokes_dv_boundary_fptr)(const Real , const Real , const Real);
  bool (*stokes_is_dv_boundary_fptr)(const Real , const Real , const Real);
  Real (*stokes_dp_boundary_fptr)(const Real , const Real , const Real);
  bool (*stokes_is_dp_boundary_fptr)(const Real , const Real , const Real);
  Number (*initial_lvlset_fptr)(const Point&, const Parameters&, const std::string&,const std::string&);


  //AutoPtr<LinearImplicitSystem> system_density;
  //LinearImplicitSystem & system2;

  

public:

    ~Multiphase() {
        delete exodusobj;
    };
    Multiphase(int argc, char** argv);//,Mesh& meshref);
    bool simulating();
    void set_n_timesteps(unsigned int numberOfTimesteps);
    int get_n_timesteps();
    void initialize(GetPot& input_file);
    void project_solutions();
    void stop();
    void increase_timestep();
    void initial_refinement();
    bool refinement(int ref_steps,bool refine_by_stokes_error=1,bool refine_interface=0);
    bool refinement(void);
    void level_reinitialization(int numberOfSteps);
    void level_reinitialization(void);
    void calculate_curvature_heaviside();
    void solve_navier_stokes();
    ErrorVector calculate_stokes_error ();
    void level_advection();
    void save_global_results();
    void save_results_to_file();
    void calculate_timestep();
    void print_cfl();
    void store_ns_old_solutions();
    void attach_initial_lvlset_function(Number fptr(const Point&, const Parameters&,
            const std::string&,const std::string&))
    {
             libmesh_assert (fptr != NULL);
             initial_lvlset_fptr = fptr;
    }

    void attach_stokes_v_dirichlet_functions(bool ns_is_dv_boundary_fptr(const Real,
                      const Real,
                      const Real), RealVectorValue ns_dv_boundary_fptr(const Real,
                      const Real ,
                      const Real )){
        libmesh_assert (ns_is_dv_boundary_fptr!=NULL );
        stokes_is_dv_boundary_fptr = ns_is_dv_boundary_fptr;

        libmesh_assert (ns_dv_boundary_fptr!=NULL);
        stokes_dv_boundary_fptr = ns_dv_boundary_fptr;
    }
    void attach_stokes_p_dirichlet_functions(bool ns_is_dp_boundary_fptr(const Real,
                      const Real,
                      const Real), Real ns_dp_boundary_fptr(const Real,
                      const Real ,
                      const Real )){
        libmesh_assert (ns_is_dp_boundary_fptr!=NULL );
        stokes_is_dp_boundary_fptr = ns_is_dp_boundary_fptr;

        libmesh_assert (ns_dp_boundary_fptr!=NULL);
        stokes_dp_boundary_fptr = ns_dp_boundary_fptr;
    }



private:
    Number get_value_on_coordinate(Real x, Real y, Real z,const System& system, const std::string& unknown_name);
    Number initial_stokes_function (const Point& p,
                      const Parameters& parameters,
                      const std::string& sys_name,
                      const std::string& var_name)
    {
      RealVectorValue vec = dirichlet_boundary_v_stokes2(p(0),p(1),p(2));

      if(var_name=="u")
      {
          return vec(0);
      }
      if(var_name=="v")
      {
          return vec(1);
      }
      if(var_name=="w")
      {
          return vec(1);
      }
    return 0;
  }

    

};

