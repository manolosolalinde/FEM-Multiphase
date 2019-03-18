//Aseembles for ex3_2

//#ifndef _ASSEMBLE_SYSTEMS_H
//#define	_ASSEMBLE_SYSTEMS_H

#include "iostream"
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"
#include "exact_solution.h"
#include "mesh_data.h"
#include "perf_log.h"
#include "boundary_info.h"
#include "utility.h"


// Define the Finite Element object.
#include "fe.h"

// Define Gauss quadrature rules.
#include "quadrature_gauss.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"


// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

#include "elem.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"
#include "transient_system.h"
#include <exodusII_io.h>
#include "interface.h"
#include "quadrature_rules.h"
#include "quadrature.h"
#include "nonlinear_solver.h"
#include "nonlinear_implicit_system.h"

#include "advdiffsystem.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"




void assemble_kurvature (EquationSystems& es,
                  const std::string& system_name);

void assemble_aux_stokes (EquationSystems& es,
                      const std::string& system_name);

void assemble_stokes (EquationSystems& es,
                      const std::string& system_name);

void assemble_normal (EquationSystems& es,
                      const std::string& system_name);
void assemble_heaviside (EquationSystems& es,
                      const std::string& system_name);

void assemble_corrector (EquationSystems& es,
                  const std::string& system_name);

void correct_reinit (EquationSystems& es);

void assemble_refinement (EquationSystems& es,
                  const std::string& system_name);

void assemble_corrector (EquationSystems& es,
                  const std::string& system_name);
void assemble_test (EquationSystems& es,
                  const std::string& system_name);


RealVectorValue normalx_velocity (EquationSystems& es,
              const std::string& system_name);

Real normalx_force (EquationSystems& es,
              const std::string& system_name);

RealVectorValue normaly_velocity (EquationSystems& es,
              const std::string& system_name);

Real normaly_force (EquationSystems& es,
              const std::string& system_name);

RealVectorValue reinit_velocity (EquationSystems& es,
              const std::string& system_name);

Real reinit_force (EquationSystems& es,
              const std::string& system_name);

Real reinit_epsilon (EquationSystems& es,
              const std::string& system_name);

RealVectorValue lvlset_velocity (EquationSystems& es,
              const std::string& system_name);

void compute_jacobian (const NumericVector<Number>& soln,
                       SparseMatrix<Number>&  jacobian,
                       NonlinearImplicitSystem& sys);

void compute_residual (const NumericVector<Number>& soln,
                       NumericVector<Number>& residual,
                       NonlinearImplicitSystem& sys);

void assemble_jump (EquationSystems& es,
                      const std::string& system_name);

void assemble_force (EquationSystems& es,
                      const std::string& system_name);

void assemble_whole_force (EquationSystems& es,
                      const std::string& system_name);

void compute_navier_stokes (const NumericVector<Number>& soln,
                       SparseMatrix<Number>&  jacobian,
                       NumericVector<Number>& residual,
                       NonlinearImplicitSystem& sys,
                       unsigned int flag);

