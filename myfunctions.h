// C++ Includes
#include <math.h>

// Mesh library includes
#include "libmesh_common.h"
#include "point.h"
#include "parameters.h"
#include "vector_value.h"
#include "utility.h"


Number function_eliptic_bubble (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&);

Number function_one_bubble (const Point& p,
                const Parameters& parameters,
                const std::string&,
                const std::string&);

Number function_two_bubbles (const Point& p,
                const Parameters& parameters,
                const std::string&,
                const std::string&);
Number function_testreinit (const Point& p,
                const Parameters& parameters,
                const std::string&,
                const std::string&);

Number function_supg (const Point& p,
                const Parameters& parameters,
                const std::string&,
                const std::string&);
Number function_one_phase (const Point& p,
                const Parameters& parameters,
                const std::string&,
                const std::string&);
Number exact_0_solution(const Point& p,
                  const Parameters&,   // EquationSystem parameters, not needed
                  const std::string&,  // sys_name, not needed
                  const std::string&); // unk_name, not needed);

Gradient exact_0_derivative(const Point& p,
                      const Parameters&,   // EquationSystems parameters, not needed
                      const std::string&,  // sys_name, not needed
                      const std::string&); // unk_name, not needed);

RealVectorValue dirichlet_boundary0(const Real x,
                      const Real y,
                      const Real z=0.);

bool is_dirichlet_boundary_v_stokes1(const Real x,
                      const Real y,
                      const Real z=0.);

RealVectorValue dirichlet_boundary_v_stokes1(const Real x,
                      const Real y,
                      const Real z=0.);

bool is_dirichlet_boundary_v_stokes2(const Real x,
                      const Real y,
                      const Real z=0.);

RealVectorValue dirichlet_boundary_v_stokes2(const Real x,
                      const Real y,
                      const Real z=0.);

bool is_dirichlet_boundary_p_stokes1(const Real x,
                      const Real y,
                      const Real z);

Real dirichlet_boundary_p_stokes1(const Real x,
                      const Real y,
                      const Real z=0.);


Real dirichlet_boundary2(const Real x,
                      const Real y,
                      const Real z=0.);

bool is_dirichlet_boundary_all(const Real x,
                      const Real y,
                      const Real z=0.);

bool is_dirichlet_boundary_1point(const Real x,
                      const Real y,
                      const Real z=0.);


