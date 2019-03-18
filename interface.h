/* 
 * File:   Interface.h
 * Author: manolo
 *
 * Created on February 27, 2010, 10:11 PM
 */

#ifndef _INTERFACE_H
#define	_INTERFACE_H

//C++ Includes
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
#include <point_locator_tree.h>
#include "mesh_data.h"

// Define the Finite Element object.
#include "fe.h"

// Define Gauss quadrature rules.
#include "quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "elem.h"

typedef std::map< const Node*, Point >  Interface_Edge;
typedef std::map<const Elem*, std::vector<bool>   > Interface_Elem;


//
// Define the DofMap, which handles degree of freedom
// indexing.

class Interface {
    const MeshBase& _mesh;
    const EquationSystems& _equation_systems;
    const std::string _sys_name;
    Real _epsilon;
    Real _deltax;
    Real _deltax_min;
    Real _epsilon_min;
    Real _epsilon_ref;
    int _refinements;

public:
    //const MeshBase& _mesh;
    ~Interface () {};
    Interface (EquationSystems& es,
                      const std::string& system_name);
    Interface_Elem interface_elem_data;
    Interface_Edge interface_edge_data;
    Interface_Edge unique_interface_edge_data;
    void locate_points ();
    Real H (Number lvlset_value);
    Real H (Number lvlset_value, Real epsilon);
    Real dH (Real lvlset_value);
    Real dH (Real lvlset_value, Real epsilon);
    Real dH_min (Real lvlset_value);
    Real Cone_ref (Real lvlset_value, Real ceil);
    Real S (Real lvlset0_value,Real custom_delta);
    void print_info(std::ostream &out=libMesh::out);
    void set_refinement_level(unsigned int base_max_h_level);
    short int sgn(Number z0);

    Real getEpsilon_ref() const {
        return _epsilon_ref;
    }

    void setEpsilon_ref(Real _epsilon_ref) {
        this->_epsilon_ref = _epsilon_ref;
    }

    Real getEpsilon() const {
        return _epsilon;
    }

    void setEpsilon(Real _epsilon) {
        this->_epsilon = _epsilon;
    }

    Real getDeltax_min() const {
        return _deltax_min;
    }

    void setDeltax_min(Real _deltax_min) {
        this->_deltax_min = _deltax_min;
    }

    Real getEpsilon_min() const {
        return _epsilon_min;
    }

    void setEpsilon_min(Real _epsilon_min) {
        this->_epsilon_min = _epsilon_min;
    }




protected:
    struct str_vint_compare
    {
        bool operator() (const void * v1, const void * v2) {
        int *f1, *f2;
        f1 = (int*) v1;
        f2 = (int*) v2;
        return (f1[0]<f2[0]);
        }
    };
    struct str_vint_compare2
    {
        bool operator() (const void * v1, const void * v2) {
        int *f1, *f2;
        f1 = (int*) v1;
        f2 = (int*) v2;
        return (f1[0]==f2[0]);
        }
    };
    struct str_vint_compare3
    {
        bool operator() (const void * v1, const void * v2) {
        int *f1, *f2;
        f1 = (int*) v1;
        f2 = (int*) v2;
        return (f1[1]<f2[1]);
        }
    };

    
};

/*class InterfaceHelper {
public:
    static
};*/





#endif	/* _INTERFACE_H */

