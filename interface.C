/* 
 * File:   Interface.c
 * Author: manolo
 *
 * Created on February 27, 2010, 10:17 PM
 */

#include <stdio.h>
#include <stdlib.h>
//C++ Includes
#include "iostream"
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "linear_implicit_system.h"

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


#include "interface.h"


Interface::Interface(EquationSystems& es, const std::string& system_name) :
     _mesh                    (es.get_mesh()),
     _equation_systems        (es),
     _sys_name                (system_name),
     _epsilon                 (0),
     _epsilon_min             (0),
     _epsilon_ref             (0),
     _refinements             (0),
     _deltax                  (0),
     _deltax_min              (0)

{
    if(_mesh.elements_begin()!=_mesh.elements_end())
        this->set_refinement_level(1);
}

void Interface::set_refinement_level(unsigned int base_max_h_level)
{
    if (base_max_h_level == 0 )
    {
        libMesh::err << "ERROR: min base_max_h_level is 1" << std::endl;
                libmesh_error();
    }
    // OBS: base_max_h_level = 1 means no refinement
    unsigned int refinement = base_max_h_level-1;
        _refinements = refinement;
        const Elem* parent = (*(_mesh.elements_begin()))->top_parent();
        _deltax = parent->hmin();
        _deltax_min = _deltax/pow(2,_refinements);
}

//
void Interface::locate_points()
{

//initialization

const System& system4 = _equation_systems.get_system (_sys_name.c_str());

NumericVector<Number>& lvlset_vector = *system4.current_local_solution;

interface_edge_data.clear();
interface_elem_data.clear();
unique_interface_edge_data.clear();



std::vector <int*> edge2d(_mesh.n_elem()*3);
std::vector <int*> edge2d_init(_mesh.n_elem()*3);

  int ii=0;


  //MeshData mesh_data(_mesh); unnecessary
  //mesh_data.activate(); unnecessary

  MeshBase::const_element_iterator       el     = _mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.elements_end();
  // Loop over the elements.  Note that  ++el is preferred to
  // el++ since the latter requires an unnecessary temporary
  // object.
  for ( ; el != end_el ; ++el)
  {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      for (int i=3 ; i<6; i++)
      {
          edge2d_init[ii] = new int[3];
          edge2d_init[ii][0]=elem->get_node(i)->id();
          edge2d_init[ii][1]=elem->id();
          edge2d_init[ii][2]=i;
          ii++;
      }
   }// done iterating elements

  // unnecesary mesh_data.insert_elem_data(artificial_data); //store artificialdata into meshdata

  str_vint_compare vint_compare;
  str_vint_compare2 vint_compare2;
  str_vint_compare3 vint_compare3;

 
  edge2d = edge2d_init;

  std::sort(edge2d.begin(), edge2d.end(), vint_compare);

  std::vector<int*>::iterator it;
  it =  std::unique(edge2d.begin(), edge2d.end(), vint_compare2);
  edge2d.resize(it - edge2d.begin());

  /* //PRINT EDGES
  std::vector <int*> edge2d_edgetoel(_mesh.n_elem()*3);
  std::vector <int*> edge2d_eltoedge(_mesh.n_elem()*3);

  edge2d_eltoedge = edge2d;
  edge2d_edgetoel = edge2d;
  std::sort(edge2d_eltoedge.begin(), edge2d_eltoedge.end(), vint_compare3);
  for (int i=0;i<25;i++)
  {
      //std::cout << edge2d_eltoedge[i][0] << "   " << edge2d_eltoedge[i][1]<< "   "<< edge2d_eltoedge[i][2] <<std::endl;
      std::cout << edge2d_edgetoel[i][0] << "   " << edge2d_edgetoel[i][1]<< "   "<< edge2d_edgetoel[i][2] <<std::endl;
  }*/


  //Add interface points algorithm

  
  std::vector<Point*> interface_points (edge2d.size()*2); //stores interface points
  std::vector<int*> inter_edge (interface_points.size());//edge ids for the interface points
  int count = 0; // interface points counter;
  const Node *za, *zb, *zc;
  Point ap,ipoint;


  for (int i=0;i<edge2d.size();i++) //START ITERATING EDGES
  {

      za = _mesh.node_ptr(edge2d[i][0]);
      
      if (za->id()!=edge2d[i][0])
          std::cout << "ESTO NO PUEDE SUCEDER. za->id()!=edge2d[i][0]";

      //assigning pointers to nodes of the edge
      if (edge2d[i][2]==3)
      {
          zb = _mesh.elem(edge2d[i][1])->get_node(0);
          zc = _mesh.elem(edge2d[i][1])->get_node(1);
      }
      else
      {
          if (edge2d[i][2]==4)
          {
              zb = _mesh.elem(edge2d[i][1])->get_node(1);
              zc = _mesh.elem(edge2d[i][1])->get_node(2);
          }
          else
          {
          if (edge2d[i][2]==5)
          {
              zb = _mesh.elem(edge2d[i][1])->get_node(2);
              zc = _mesh.elem(edge2d[i][1])->get_node(0);
          }
          else
          {
              libmesh_error();
              //manolo error;
          }
          }
      }//done assigning pointers to nodes of the edge

      short int sgnza,sgnzb,sgnzc;

      Number za_value = lvlset_vector(za->id());
      Number zb_value = lvlset_vector(zb->id());
      Number zc_value = lvlset_vector(zc->id());
      sgnza = sgn(za_value);
      sgnzb = sgn(zb_value);
      sgnzc = sgn(zc_value);

      //to be modified!!!
      if (std::abs(sgnza+sgnzb)<2) //we've got interface between za and zb
      {
          if (sgnzb ==0)
          {
              ap = *zb - *za;
          }
          else
          {
              if ((sgnza+sgnzb==0) && (sgnza!=0))
              {
                  ap = za_value/(za_value-zb_value)*(*zb-*za);
              }
          }
          if (sgnza==0)
          {
             std::cout << "I was hwere";
             ap.zero();
          }

        interface_points[count] = new Point;
        interface_points[count][0] = *za + ap; //include interface point
        inter_edge[count] = new int[3];
        inter_edge[count][0] = za->id();
        inter_edge[count][1] = edge2d[i][1];//element
        inter_edge[count][2] = edge2d[i][2]-3;//element side corresponding to edge
        count++;

      }
      else
      {
        if (std::abs(sgnza+sgnzc)<2) //we've got interface between za and zb
         {
          if (sgnzc ==0)
          {
              ap = *zc - *za;
          }
          else
          {
              if ((sgnza+sgnzc==0) && (sgnza!=0))
              {
                  ap = za_value/(za_value-zc_value)*(*zc-*za);
              }
          }
          if (sgnza==0)
          {
             std::cout << "I was hwere, this can never happen";
             ap.zero();
          }

        interface_points[count] = new Point;
        interface_points[count][0] = *za + ap; //include interface point
        inter_edge[count] = new int[3];
        inter_edge[count][0] = za->id();
        inter_edge[count][1] = edge2d[i][1];//element id
        inter_edge[count][2] = edge2d[i][2]-3;//element side
        count++;
       }
      }

  }//done iterating edges

  interface_points.resize(count);
  inter_edge.resize(count);


  //std::map<const Elem*, std::vector<Number> > interface_data;

  //std::map<const Elem*, std::vector<bool> > mesh_interface_data;
  //mesh_interface_data.clear();

  //std::map<const Elem*, std::vector<bool>   > interface_elem_data;
  //std::map<const Node*, Point* > interface_edge_data;

  //store the information about the interface points in each node
  for (int i=0;i<interface_points.size();i++)
  {
      int currentedgeid = inter_edge[i][0];
      const Node * currentedge = _mesh.node_ptr(currentedgeid);

      //filling interface edge point information
      interface_edge_data.insert(std::make_pair(currentedge,*interface_points[i]));
  }



  //check which elements have interface and store the information
  for (int i=0;i<interface_points.size();i++)
  {
      Elem * elema = _mesh.elem(inter_edge[i][1]);
      Elem * elemb = elema->neighbor(inter_edge[i][2]);

      //filling interface element elema information
      if (interface_elem_data.find(elema)==interface_elem_data.end())//couldnt find elema in interface_edge_daga
      {
        std::vector<bool> eleminfo(3);
        for(int i=3;i<6;i++)
        {
          const Node * currentedge = _mesh.node_ptr(elema->node(i));
          //elema
          if (interface_edge_data.find(currentedge)!=interface_edge_data.end()) //SIDE WITH INTERFACE
              eleminfo[i-3]=1;
        }
        interface_elem_data.insert(std::make_pair(elema,eleminfo));
      }

      //filling interface element elemb information
      if (elemb!=NULL)
      {
      if (interface_elem_data.find(elemb)==interface_elem_data.end())//couldnt find elema in interface_edge_daga
      {
        std::vector<bool> eleminfo(3);
        for(int i=3;i<6;i++)
        {
          const Node * currentedge = _mesh.node_ptr(elemb->node(i));
          //elema
          if (interface_edge_data.find(currentedge)!=interface_edge_data.end()) //SIDE WITH INTERFACE
              eleminfo[i-3]=1;
        }
        interface_elem_data.insert(std::make_pair(elemb,eleminfo));
      }
      }
  }

  /*unique_interface_edge_data = interface_edge_data;
  {
      Interface_Edge::iterator it3,it4;
      for(it3=unique_interface_edge_data.begin();it3!=unique_interface_edge_data.end();it3++)
      {
          for (it4=it3;it4!=unique_interface_edge_data.end();it4++)
          {
              if (it3->second==it4->second && it3!=it4)
              {
                  unique_interface_edge_data.erase(it4);
                  it4--;
              }
          }
      }
  }*/



  /*std::cout << "PRINTING INTERFACE POINTS: "<< count << std::endl;
  for (int i=0;i<count;i++)
  {
      std::cout << "EDGE NUMBER: " << inter_edge[i][0] <<"  ";
      interface_points[i][0].print(std::cout);
  }*/
  {
  /*std::cout << "PRINTING INTERFACE EDGES: "<< interface_edge_data.size() << std::endl;
  Interface_Edge::iterator it2 = interface_edge_data.begin();
  const Interface_Edge::iterator end_it = interface_edge_data.end();

  for ( ; it2!=end_it;++it2)
  {
      (it2->first)->print(std::cout);
  }*/
  }

  /*{
  std::cout << "PRINTING UNIQUE INTERFACE EDGES: "<< unique_interface_edge_data.size() << std::endl;
  Interface_Edge::iterator it2 = unique_interface_edge_data.begin();
  const Interface_Edge::iterator end_it = unique_interface_edge_data.end();

  for ( ; it2!=end_it;++it2)
  {
      (it2->first)->print(std::cout);
  }
  }*/

  // LOCATE INTERFACE POINTS RESULTS:
    // interface_elem_data: map std::map<Elem*, std::vector<bool>(3)> the vector contains info of weather the side contains interface points or not
    // interface_edge_data: map std::map<const Node*, Point* > contains pointers to nodes that have interface w/ interface location

  //AVOIDING MEMORY LEAKS
  for (int i=0;i<edge2d_init.size();i++)
  {
      delete[] edge2d_init[i];
  }
  for (int i=0;i<inter_edge.size();i++)
  {
      delete[] inter_edge[i];
  }
  for (int i=0;i<interface_points.size();i++)
  {
      delete interface_points[i];
  }


}

//Heaviside Function for the Level Set equation
Real Interface::H (Number lvlset_value)
{
    Real output;

    if (lvlset_value <= -_epsilon)
        return output = 0;

    if (lvlset_value >= -_epsilon && lvlset_value <= +_epsilon)
        return output = 0.5 + lvlset_value/2/_epsilon + 1/(2*libMesh::pi)*sin(libMesh::pi*lvlset_value/_epsilon);

    if (lvlset_value >= _epsilon)
        return output = 1;

    if (lvlset_value == +INFINITY )
        return output = 1;

    if (lvlset_value == -INFINITY)
        return output = 0;

    //Die a terrible death here
    return lvlset_value;
    libmesh_error();
}
Real Interface::H (Number lvlset_value, Real epsilon)
{
    Real output;

    if (lvlset_value < -epsilon)
        return output = 0;

    if (lvlset_value >= -epsilon && lvlset_value <= +epsilon)
        return output = 0.5 + lvlset_value/2/epsilon + 1/(2*libMesh::pi)*sin(libMesh::pi*lvlset_value/epsilon);

    if (lvlset_value > epsilon)
        return output = 1;

    if (lvlset_value == +INFINITY )
        return output = 1;

    if (lvlset_value == -INFINITY)
        return output = 0;

    //Die a terrible death here
    return lvlset_value;
    libmesh_error();
    

}

//Derivative dH/d(phi) of the Heaviside Function for the Level Set equation
Real Interface::dH (Real lvlset_value)
{
    Real output;

    if (lvlset_value < -_epsilon)
        return output = 0;

    if (lvlset_value >= -_epsilon && lvlset_value <= +_epsilon)
        return output = .5/_epsilon + 0.5/_epsilon*cos(libMesh::pi*lvlset_value/_epsilon);

    if (lvlset_value > _epsilon)
        return output = 0;

    return 0;
    //Die a terrible death here
   // libmesh_error();
}

Real Interface::dH (Real lvlset_value, Real epsilon)
{
    Real output;

    if (lvlset_value < -epsilon)
        return output = 0;

    if (lvlset_value >= -epsilon && lvlset_value <= +epsilon)
        return output = .5/epsilon + 0.5/epsilon*cos(libMesh::pi*lvlset_value/epsilon);

    if (lvlset_value > epsilon)
        return output = 0;

    return 0;
    //Die a terrible death here
   // libmesh_error();
}

Real Interface::Cone_ref (Real lvlset_value,Real ceil)
{
    Real output;

    if (lvlset_value < -_epsilon_ref)
        return output = 0;

    if (lvlset_value >= -0 && lvlset_value <= +_epsilon_ref)
        return output = (_epsilon_ref - lvlset_value)*ceil/_epsilon_ref;

    if (lvlset_value >= -_epsilon_ref && lvlset_value <= +0)
        return output = (_epsilon_ref + lvlset_value)*ceil/_epsilon_ref;

    if (lvlset_value > _epsilon_ref)
        return output = 0;

    return 0;
    //Die a terrible death here
   // libmesh_error();
}

Real Interface::dH_min (Real lvlset_value)
{
    Real output;
    Real alpha = 0/(2*_epsilon_min);
    Real beta = 1;

    if (lvlset_value < -_epsilon_min)
        return output = 0;

    if (lvlset_value >= -_epsilon_min && lvlset_value <= +_epsilon_min)
        return output = alpha + beta*(.5/_epsilon_min + 0.5/_epsilon_min*cos(libMesh::pi*lvlset_value/_epsilon_min));

    if (lvlset_value > _epsilon_min)
        return output = 0;

    return 0;
    //Die a terrible death here
   // libmesh_error();
}

//Recomended sustitution for sgn()
Real Interface::S (Real lvlset0_value,Real custom_delta)
{
    Real output = lvlset0_value/sqrt(pow(lvlset0_value,2)+pow(custom_delta,2));
    return output;
}


short int Interface::sgn (Number z0)
        {
        short int var;
            if (z0>0)
        return var = 1;
         if (z0<0)
        return var = -1;
          if (z0==0)
        return var = 0;
        }

void Interface::print_info(std::ostream &out){
    out << "Interface class information\n";
    out << "  Reference element size: deltax = " << _deltax <<std::endl;
    out << "  Selected h refinement level: hlevel = "<< _refinements+1 <<std::endl;
    out << "  Minimum reference element size: deltax_min = " << _deltax_min <<std::endl;
    out << "  Standand interface jump smearing: epsilon = " << _epsilon <<std::endl;
    out << "  Minimum interface jump smearing: epsilon_min = " << _epsilon_min <<std::endl;
    out << "  Refinement interface jump smearing: epsilon_ref = "<< _epsilon_ref<<std::endl;
    out <<std::endl;
}