
#ifndef __advdiffsystem_h__
#define __advdiffsystem_h__


#include "linear_implicit_system.h"
#include "numeric_vector.h"
#include "transient_system.h"
#include <auto_ptr.h>
#include "fem_system.h"

#include "quadrature_rules.h"
#include "quadrature.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "elem.h"
#include "fe.h"
#include "utility.h"
#include "equation_systems.h"
#include "libmesh.h"
#include "mesh.h"


using namespace libMesh;


struct AdvDiffStruct {
          std::string struct_name;
          RealVectorValue velocity;
          Real force;
          Real epsilon;
          Real theta;
          Real dt;
          Real time;
};


//Forward definitions
//class libMesh::TransientSystem;
//class libMesh::LinearImplicitSystem;

template <class Base>
class AdvectionDiffusionSystem : public Base
{
public:
    AdvectionDiffusionSystem(EquationSystems& es,
               const std::string& name,
               const unsigned int number);

    virtual ~AdvectionDiffusionSystem();

    virtual void clear ();

    virtual void reinit ();

    Elem* elem_ptr;

    AdvDiffStruct data;

    //data_var_list_all
    std::map< std::string, std::vector<unsigned int> > data_var_list_all;
    
    //data_systems Stores the system numbers
    std::vector<unsigned int> data_systems;

    //data_var_list
    //lista de las variables que no pertenecen al systema
    std::vector<std::string> data_var_list;

    std::map<std::string,Real*> data_var_map;
    std::map<std::string,Gradient*> data_var_map_gradient;
    std::map<std::string,Real> data_var_map_mean;
    std::map<std::string,Gradient> data_var_map_mean_gradient;
    RealVectorValue (* f_velocity_ptr) (EquationSystems& es,
                  const std::string& system_name);
    Real (* f_force_ptr) (EquationSystems& es,
                  const std::string& system_name);
    Real (* f_epsilon_ptr) (EquationSystems& es,
                  const std::string& system_name);
    
    void use_on_assemble(System& sys);
    void attach_variable_for_assemble(System& sys,const std::string& var_name,const std::string& grad="base");
    void attach_velocity_function(RealVectorValue fptr(EquationSystems& es,
                  const std::string& sys_name));
    void attach_force_function(Real fptr(EquationSystems& es,
                  const std::string& sys_name));
    void attach_epsilon_function(Real fptr(EquationSystems& es,
                  const std::string& sys_name));

    Number get_variable (std::string var_name);

    Gradient get_gradient (std::string grad_name);

    Elem* get_element_ptr ();

    typedef AdvectionDiffusionSystem<Base> sys_type;

    sys_type & system () { return *this; }

    Number old_solution (const unsigned int global_dof_number) const;


    virtual std::string system_type() const;

    static Real bubble1(const Point& qpoint){

        Real x = qpoint(0);
        Real y = qpoint(1);
        //TODO this only works with QUAD elements
        if (x<=1 && x>=-1 && y<=1 && y>=-1)
        {
            return (1-x*x)*(1-y*y);
        }
        else
        {
            libmesh_error();
        }

        libmesh_error();
    }

    static Gradient grad_bubble1(const Point& qpoint){

        Gradient mygrad;
        Real x = qpoint(0);
        Real y = qpoint(1);
        //TODO this only works with QUAD elements
        if (x<=1 && x>=-1 && y<=1 && y>=-1)
        {
            mygrad(0)=2*x*(y*y-1);
            mygrad(1)=2*(x*x-1)*y;
            return mygrad;
        }
        else
        {
            libmesh_error();
        }

        libmesh_error();

    }

    static Real bubble2(const Point& qpoint){

        Real x = qpoint(0);
        Real y = qpoint(1);
        int flag = 0;//region flag
        const Real Xb = -0.99;
        const Real Yb = -0.99;

        //TODO this only works with QUAD elements
        if (x<=1 && x>=-1 && y<=1 && y>=-1)
        {
            if (x>=Xb && y>=Yb)//Region 3 or 2
            {
                if (y >= Yb+ (x-Xb)*(1-Yb)/(1-Xb)) //region 3
                    flag = 3;
                else
                    flag = 2;
            }
            if (x<=Xb && y>=Yb)//Region 3 or 4
            {
                if (y>=Yb + (Xb-x)*(1-Yb)/(Xb+1)) //region 3
                    flag = 3;
                else
                    flag = 4;
            }
            if (x<=Xb && y<=Yb)//Region 4 or 1
            {
                if (y>= Yb - (Xb-x)*(Yb+1)/(Xb+1))//Region 4
                    flag = 4;
                else
                    flag = 1;
            }
            if (x>=Xb && y<=Yb)//Region 1 or 2
            {
                if (y>=Yb-(x-Xb)*(Yb+1)/(1-Xb))
                    flag = 2;
                else
                    flag =1;
            }

            switch (flag)
            {
                case 1:
                    return x/(Xb*(1+Yb))*(1+y);
                    break;
                case 2:
                    return y/(Yb*(1-Xb))*(1-x);
                    break;
                case 3:
                    return x/(Xb*(1-Yb))*(1-y);
                    break;
                case 4:
                    return y/(Yb*(1+Xb))*(1+x);
                    break;
                default:
                    libmesh_error();//none of the cases matched
            }


        }
        else
        {
            libmesh_error();
        }

        libmesh_error();
    }

    static Gradient grad_bubble2(const Point& qpoint){

        Gradient mygrad;
        Real x = qpoint(0);
        Real y = qpoint(1);
        int flag = 0;//region flag
        const Real Xb = -0.99;
        const Real Yb = -0.99;

        //TODO this only works with QUAD elements
        if (x<=1 && x>=-1 && y<=1 && y>=-1)
        {
            if (x>=Xb && y>=Yb)//Region 3 or 2
            {
                if (y >= Yb+ (x-Xb)*(1-Yb)/(1-Xb)) //region 3
                    flag = 3;
                else
                    flag = 2;
            }
            if (x<=Xb && y>=Yb)//Region 3 or 4
            {
                if (y>=Yb + (Xb-x)*(1-Yb)/(Xb+1)) //region 3
                    flag = 3;
                else
                    flag = 4;
            }
            if (x<=Xb && y<=Yb)//Region 4 or 1
            {
                if (y>= Yb - (Xb-x)*(Yb+1)/(Xb+1))//Region 4
                    flag = 4;
                else
                    flag = 1;
            }
            if (x>=Xb && y<=Yb)//Region 1 or 2
            {
                if (y>=Yb-(x-Xb)*(Yb+1)/(1-Xb))
                    flag = 2;
                else
                    flag =1;
            }

            switch (flag)
            {
                case 1:
                    mygrad(0) = (y+1)/(Xb*Yb+Xb);
                    mygrad(1) = x/(Xb*Yb+Xb);
                    return mygrad;
                    break;
                case 2:
                    mygrad(0) = y/((Xb-1)*Yb);
                    mygrad(1) = (x-1)/((Xb-1)*Yb);
                    return mygrad;
                    break;
                case 3:
                    mygrad(0) = (y-1)/(Xb*(Yb-1));
                    mygrad(1) = x/(Xb*(Yb-1));
                    return mygrad;
                    break;
                case 4:
                    mygrad(0) = y/(Xb*Yb+Yb);
                    mygrad(1) = (x+1)/(Xb*Yb+Yb);
                    return mygrad;
                    break;
                default:
                    libmesh_error();//none of the cases matched
            }


        }
        else
        {
            libmesh_error();
        }

        libmesh_error();
    }




    


private:

    EquationSystems & _es;

    static void assemble_advdiff (EquationSystems& es,
                  const std::string& system_name);

    static System& systems (AdvectionDiffusionSystem& mysys, unsigned int i);
    static DofMap& dof_map_systems (AdvectionDiffusionSystem& mysys, unsigned int i);

    AdvectionDiffusionSystem* sys_ptr;


    std::vector<std::string> data_var_list_mean;

    std::vector<std::string> data_var_list_mean_gradient;

    std::vector<bool> data_myvar_check;
    //GUIDE FOR data_myvar_check
        //data_myvar_check[0]  <->  mysystem.oldsolution
        //data_myvar_check[1]  <->  gradient(mysystem.oldsolution)
        //data_myvar_check[2]  <->  mean(mysystem.oldsolution)
        //data_myvar_check[3]  <->  mean(gradient(mysystem.oldsolution))
        //data_myvar_check[4]  <->  d2dx2(mysystem.oldsolution)
        //data_myvar_check[5]  <->  d2dy2(mysystem.oldsolution)


      
};

typedef AdvectionDiffusionSystem<TransientLinearImplicitSystem> TransientAdvectionDiffusionSystem;
typedef AdvectionDiffusionSystem<LinearImplicitSystem> StaticAdvectionDiffusionSystem;


// AdvectionDiffusionSystem inline methods
template <class Base>
inline
std::string AdvectionDiffusionSystem<Base>::system_type () const
{
std::string type = "AdvectionDiffusion";
type += Base::system_type ();

return type;
}

#endif