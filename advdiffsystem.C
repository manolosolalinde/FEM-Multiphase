#include "advdiffsystem.h"

#include "assemble_systems.h"




//Class constructor
template <class Base>
AdvectionDiffusionSystem<Base>::AdvectionDiffusionSystem(EquationSystems& es,
               const std::string& name,
               const unsigned int number) : Base(es, name, number),_es(es)
    {
        this->attach_assemble_function(assemble_advdiff);
        data.struct_name = name;
        data.struct_name += " DATA";
        data.epsilon = 0;
        data.theta = 1;
        data.dt = 0;
        data.time = -10;
        data.velocity = NumberVectorValue(0,0);//FIXME - Not ready for 3d
        data.force = 0;
        es.parameters.set<AdvDiffStruct*>(data.struct_name) = &data;
        sys_ptr = this;
        data_var_list.clear();
        data_var_list_all.clear();
        data_systems.clear();
        f_epsilon_ptr = NULL;
        f_force_ptr = NULL;
        f_velocity_ptr = NULL;
        elem_ptr = NULL;
        data_myvar_check.resize(6,0);
        data_myvar_check[0]=1;
        data_myvar_check[1]=1;
   
    }

template <class Base>
AdvectionDiffusionSystem<Base>::~AdvectionDiffusionSystem()
{
    this->clear();
}

template <class Base>
void AdvectionDiffusionSystem<Base>::clear ()
{
    // clear the parent data
    Base::clear();
}

template <class Base>
void AdvectionDiffusionSystem<Base>::reinit ()
{
    Base::reinit();
}



template <class Base>
void AdvectionDiffusionSystem<Base>::attach_variable_for_assemble(System& sys,const std::string& var_name,const std::string& grad)
{
//    std::vector<std::string> var_names;
//    _es->build_variable_names(var_names);

    libmesh_assert(sys.variable_number(var_name)>=0);
    
    //Note that variables on the lists below may be repeated
    unsigned int switcher=0;
    if (grad.compare("base")==0)
        switcher = 1;
    if (grad.compare("grad")==0)
        switcher = 2;
    if (grad.compare("mean")==0)
        switcher = 3;
    if (grad.compare("mean_grad")==0)
        switcher = 4;
    if (grad.compare("old")==0)
        switcher = 5;
    if (grad.compare("old_grad")==0)
        switcher = 6;
    if (grad.compare("old_mean")==0)
        switcher = 7;
    if (grad.compare("old_grad_mean")==0)
        switcher = 8;
    if (grad.compare("d2dx2")==0)
        switcher = 9;
    if (grad.compare("d2dy2")==0)
        switcher = 10;

    //GUIDE FOR data_myvar_check
        //data_myvar_check[0]  <->  mysystem.oldsolution
        //data_myvar_check[1]  <->  gradient(mysystem.oldsolution)
        //data_myvar_check[2]  <->  mean(mysystem.oldsolution)
        //data_myvar_check[3]  <->  mean(gradient(mysystem.oldsolution))
        //data_myvar_check[4]  <->  d2dx2(mysystem.oldsolution)
        //data_myvar_check[5]  <->  d2dy2(mysystem.oldsolution)


    //if the variable belongs to this system
    if (var_name.compare(this->variable_name(0))==0)
    {
        if (switcher <= 4)
        {
            libMesh::err << "ERROR: You can not attach "
                      << var_name <<", "<< grad
                      << "to the system." << std::endl;
            libmesh_error();
        }
        if (switcher == 5 || switcher == 6)
        {
            libMesh::err << "WARNING: You don't need to attach "
                      << var_name <<", "<< grad
                      << "to the system." << " Use 'uold' to call "
                      << "this variable." << std::endl;
        }

        if (switcher == 7)
            data_myvar_check[2]=1;

        if (switcher == 8)
            data_myvar_check[3]=1;

        if (switcher == 9)
            data_myvar_check[4]=1;

        if (switcher == 10)
            data_myvar_check[5]=1;
    }
    else //the variable belong to another system
    {

    std::vector<std::string>::iterator it;
    it = std::find(data_var_list.begin(),data_var_list.end(),var_name);
    if (data_var_list.end() == it)
    {
        //Now add the new system to the list
        std::vector<unsigned int>::iterator it2;
        it2 = std::find(data_systems.begin(),data_systems.end(),sys.number());
        if (data_systems.end() == it2)
            data_systems.push_back(sys.number());
        
        //Add the new variable to the list
        data_var_list.push_back(var_name);
        std::vector<unsigned int> var_info(4,0);

        var_info[0]=sys.number();
        var_info[1]=sys.variable_number(var_name);
        //data_var_list_all.insert(std::pair< std::string, std::vector<unsigned int> >(var_name,var_info));
        data_var_list_all[var_name] = var_info;
    }

    switch (switcher)
    {
        case 1:
            data_var_list_all[var_name][2]=1;
            break;
        case 2:
            data_var_list_all[var_name][3]=1;
            break;
        case 3:
            data_var_list_mean.push_back(var_name);
            break;
        case 4:
            data_var_list_mean_gradient.push_back(var_name);
            break;
        default:
            libMesh::err << "ERROR: Input string '" << grad <<"' is not valid. "
                         << "use 'base','grad','mean' or 'mean_grad' instead" <<std::endl;
            libmesh_error();
    }

    }

    //Print information
    /*std::map< std::string, std::vector<unsigned int> >::iterator iter;
    for ( iter=data_var_list_all.begin() ; iter != data_var_list_all.end(); iter++ )
    {
      std::cout << "data_var_list_all[" << (*iter).first.c_str() << "][0] = "
              << data_var_list_all[(*iter).first][0] <<std::endl;
      std::cout << "data_var_list_all[" << (*iter).first.c_str() << "][1] = "
              << data_var_list_all[(*iter).first][1] <<std::endl;
      std::cout << "data_var_list_all[" << (*iter).first.c_str() << "][2] = "
              << data_var_list_all[(*iter).first][2] <<std::endl;
      std::cout << "data_var_list_all[" << (*iter).first.c_str() << "][3] = "
              << data_var_list_all[(*iter).first][3] <<std::endl;

    }
    std::cout << "end of list " <<std::endl;*/



}

template <class Base>
void AdvectionDiffusionSystem<Base>::attach_velocity_function(RealVectorValue fptr(EquationSystems&,const std::string&))
{
     libmesh_assert (fptr != NULL);
     f_velocity_ptr = fptr;
}

template <class Base>
void AdvectionDiffusionSystem<Base>::attach_force_function(Real fptr(EquationSystems&,const std::string&))
{
     libmesh_assert (fptr != NULL);
     f_force_ptr = fptr;
}

template <class Base>
void AdvectionDiffusionSystem<Base>::attach_epsilon_function(Real fptr(EquationSystems&,const std::string&))
{
     libmesh_assert (fptr != NULL);
     f_epsilon_ptr = fptr;
}

template <class Base>
Number AdvectionDiffusionSystem<Base>::get_variable(std::string var_name)
{
    std::vector<std::string>::iterator it;
    it = std::find(data_var_list.begin(),data_var_list.end(),var_name);
    if (data_var_list.end() == it)
    {
        libMesh::err << "ERROR: variable '"<<var_name
                     <<"' was not attached to system '"
                     <<this->name()<< "'"<<std::endl;
        libmesh_error();
    }
    else
    {
        if(data_var_list_all[var_name][2]!=1)
        {
        libMesh::err << "ERROR: variable '"<<var_name
                     <<"' was not attached to system '"
                     <<this->name()<< "'"<<std::endl;
        libmesh_error();
        }
    }

    return *(data_var_map[var_name]);
}

template <class Base>
Gradient AdvectionDiffusionSystem<Base>::get_gradient(std::string grad_name)
{
    std::string var_name = grad_name;
    var_name.erase(0,5);
    std::vector<std::string>::iterator it;
    it = std::find(data_var_list.begin(),data_var_list.end(),var_name);
    if (data_var_list.end() == it)
    {
        libMesh::err << "ERROR: variable '"<<var_name
                     <<"' was not attached to system '"
                     <<this->name()<< "'"<<std::endl;
        libmesh_error();
    }
    else
    {
        if(data_var_list_all[var_name][3]!=1)
        {
        libMesh::err << "ERROR: variable gradient '"<<var_name
                     <<"' was not attached to system '"
                     <<this->name()<< "'"<<std::endl;
        libmesh_error();
        }
    }


    return *(data_var_map_gradient[grad_name]);
}

template <class Base>
Elem* AdvectionDiffusionSystem<Base>::get_element_ptr()
{
    return elem_ptr;
}

template <class Base>
System& AdvectionDiffusionSystem<Base>::systems (AdvectionDiffusionSystem<Base>& mysys,unsigned int sys_number)
{
    EquationSystems& es = mysys.get_equation_systems();
    System& sys = es.get_system(sys_number);
    return sys;
}

template <class Base>
DofMap& AdvectionDiffusionSystem<Base>::dof_map_systems(AdvectionDiffusionSystem<Base>& mysys, unsigned int sys_number)
{
    DofMap& dofmap = systems(mysys,sys_number).get_dof_map();
    return dofmap;
}

template<>
Number AdvectionDiffusionSystem<LinearImplicitSystem>::old_solution(const unsigned int global_dof_number) const
{
    return 0;
}

template<>
Number AdvectionDiffusionSystem<TransientLinearImplicitSystem>::old_solution(const unsigned int global_dof_number) const
{
    // Check the sizes
    libmesh_assert (global_dof_number < this->get_dof_map().n_dofs());
    libmesh_assert (global_dof_number < old_local_solution->size());

    return (*old_local_solution)(global_dof_number);
}

template <class Base>
void AdvectionDiffusionSystem<Base>::assemble_advdiff (EquationSystems& es,
                  const std::string& system_name)
{

  
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to each system object.

  System & sys0 = es.get_system(system_name);



  
  AdvectionDiffusionSystem<Base> & system =
            es.get_system< AdvectionDiffusionSystem<Base> >(system_name);

  AdvectionDiffusionSystem<Base> & sys = system;



//  System& (AdvectionDiffusionSystem::*systems)(unsigned int) = NULL;
//  systems = &AdvectionDiffusionSystem.systems;

  //pt2Object = (AdvectionDiffusionSystem*) system.sys_ptr;


//  std::cout << "SYSTEM TYPE = "<< system.system_type() <<std::endl;
//  std::cout << "SYSTEM TYPE = "<< sys0.system_type() <<std::endl;





  // Numeric ids corresponding to each variable in the systems
  //not here anymore

  FEType fe_type = system.variable_type(0);

  AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  QGauss qrule (dim,   fe_type.default_quadrature_order());
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  fe->attach_quadrature_rule      (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real>& JxW      = fe->get_JxW();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<Real> >& d2phidx2 = fe->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe->get_d2phidy2();
  const std::vector<std::vector<Real> >& psi = fe_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  //used in masud
  const std::vector<Point>& qpoints = qrule.get_points();

  const DofMap& dof_map = system.get_dof_map();

  // use dof_map_systems(i) to call a reference to each system

  //to be deleted - wrong sintax (ptr to ref not valid)
//  std::vector<DofMap*> dof_map_systems;
//  for (unsigned int i=0; i<system.data_systems.size(); i++)
//      dof_map_systems.push_back(&(systems(i).get_dof_map()));

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  std::map<std::string , std::vector<unsigned int> > dof_indices_name;
  //dof_indices_name.resize(system.data_var_list_all.size());

  Real dt;
  Real time;
  if (system.data.dt!=0)
      dt = system.data.dt;
  else
      dt = es.parameters.get<Real>   ("dt");

  if (system.data.time < 0)
      time = system.data.time;
  else
      time = es.parameters.get<Real> ("time");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();



  AdvDiffStruct * data = &system.data;
  RealVectorValue velocity = data->velocity;
  Real f = data->force;
  Real epsilon = data->epsilon;
  Real theta = data->theta;



  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      {
          std::map< std::string, std::vector<unsigned int> >::iterator it;
          for ( it=system.data_var_list_all.begin() ; it != system.data_var_list_all.end(); it++ )
          {
              std::string var_name = (*it).first;
              unsigned int system_number = (*it).second[0];
              unsigned int variable_number = (*it).second[1];
              const DofMap& dofmap_var = dof_map_systems(sys,system_number);
              std::vector<unsigned int> dof_indices_var;
              dofmap_var.dof_indices(elem,dof_indices_var,variable_number);
              dof_indices_name[var_name] = dof_indices_var;
          }
      }

      fe->reinit (elem);

      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());



      
      {//NAMEVAR MEANS
      std::vector<std::string>::iterator it;
      for(it=system.data_var_list_mean.begin();it!=system.data_var_list_mean.end();it++)
      {
          std::string var_name = (*it);
          Real name_mean=0;
          for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
              Real name=0;
              for (unsigned int l=0; l<phi.size(); l++)
                {
                  unsigned int sys_number = system.data_var_list_all[var_name][0];
                  name += phi[l][qp]*systems(sys,sys_number).current_solution(dof_indices_name[var_name][l]);
                }
              name_mean += JxW[qp]*name;
          }
          name_mean = name_mean/elem->volume();
          system.data_var_map_mean[var_name] = name_mean;
      }
      }


      {//NAMEVAR MEAN GRADIENTS
      std::vector<std::string>::iterator it;
      for(it=system.data_var_list_mean_gradient.begin();it!=system.data_var_list_mean_gradient.end();it++)
      {
          std::string var_name = (*it);
          Gradient grad_name_mean;
          grad_name_mean.zero();
          for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
              Gradient grad_name;
              for (unsigned int l=0; l<phi.size(); l++)
                {
                  unsigned int sys_number = system.data_var_list_all[var_name][0];
                  grad_name.add_scaled(dphi[l][qp],systems(sys,sys_number).current_solution (dof_indices_name[var_name][l]));
                }
              grad_name_mean += JxW[qp]*grad_name;
          }
          grad_name_mean = grad_name_mean/elem->volume();
          system.data_var_map_mean_gradient[var_name] = grad_name_mean;
      }
      }
      
      //MYVAR_OLD MEANS
      if (system.data_myvar_check[2])
      {
          Real myvar_mean=0;
          for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
              Real myvar=0;
              for (unsigned int l=0; l<phi.size(); l++)
              {
                  myvar += phi[l][qp]*system.old_solution  (dof_indices[l]);
              }
              myvar_mean += JxW[qp]*myvar;
          }
          myvar_mean = myvar_mean/elem->volume();
          system.data_var_map_mean["uold"]= myvar_mean;
      }

      //MYVAR_OLD GRADIENT MEANS
      if (system.data_myvar_check[3])
      {
          Gradient grad_myvar_mean;
          grad_myvar_mean.zero();
          for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
              Gradient grad_myvar;
              for (unsigned int l=0; l<phi.size(); l++)
                {
                  grad_myvar.add_scaled(dphi[l][qp],system.old_solution  (dof_indices[l]));
                }
              grad_myvar_mean += JxW[qp]*grad_myvar;
          }
          grad_myvar_mean = grad_myvar_mean/elem->volume();
          system.data_var_map_mean_gradient["uold"] = grad_myvar_mean;
      }


      //MASUD stabilization parameter calculation
      /*Real masud_int = 0;
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
          const Point& qpoint = qpoints[qp];
          masud_int += JxW[qp]*system.bubble2(qpoint);
      }*/




      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {


          // Values to hold the old solution & its gradient.
          Number   uold = 0.;
          Number velu=0;
          Number velv=0;
          Number d2uolddx2 =0;
          Number d2uolddy2 =0;
          Gradient grad_uold;
          grad_uold.zero();
          std::map<std::string , Real> name;
          std::map<std::string, Gradient> grad_name;


          for (unsigned int l=0; l<phi.size(); l++)
            {
                  if(system.data_myvar_check[0])
                      uold   += phi[l][qp]*system.old_solution  (dof_indices[l]);
                  if(system.data_myvar_check[1])
                      grad_uold.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));

                  if(system.data_myvar_check[4])
                      d2uolddx2 +=d2phidx2[l][qp]*system.old_solution (dof_indices[l]);
                  if(system.data_myvar_check[5])
                      d2uolddy2 +=d2phidy2[l][qp]*system.old_solution (dof_indices[l]);

              std::map< std::string, std::vector<unsigned int> >::iterator it;
              for ( it=system.data_var_list_all.begin() ; it != system.data_var_list_all.end(); it++ )
              {
                  std::string var_name = (*it).first;
                  if ((*it).second[2]==1)
                      name[var_name] += phi[l][qp]*systems(sys,(*it).second[0]).current_solution(dof_indices_name[var_name][l]);
                  if ((*it).second[3]==1)
                      grad_name[var_name].add_scaled(dphi[l][qp],systems(sys,(*it).second[0]).current_solution (dof_indices_name[var_name][l]));
              }

          }
          Real duolddx = grad_uold(0);
          Real duolddy = grad_uold(1);
          Real duolddz;
          if (dim==3)
              duolddz = grad_uold(2);


          //Set pointers to variables ins data_var_map
          system.data_var_map.clear();
          system.data_var_map_gradient.clear();
          std::map< std::string, std::vector<unsigned int> >::iterator it;
          for ( it=system.data_var_list_all.begin() ; it != system.data_var_list_all.end(); it++ )
          {
              std::string var_name = (*it).first;
              system.data_var_map[var_name] = &name[var_name];
              system.data_var_map_gradient["grad_" + var_name] = &grad_name[var_name];
          }

          system.data_var_map.insert(std::pair<std::string,Real*>("uold",&uold));
          system.data_var_map.insert(std::pair<std::string,Real*>("d2uolddx2",&d2uolddx2));
          system.data_var_map.insert(std::pair<std::string,Real*>("d2uolddy2",&d2uolddy2));
          system.data_var_map_gradient.insert(std::pair<std::string,Gradient*>("grad_uold",&grad_uold));

          system.elem_ptr = *el;

          if (system.f_velocity_ptr!=NULL)
          {
              velocity = system.f_velocity_ptr(es,system_name);
          }

          if (system.f_force_ptr!=NULL)
          {
              f = system.f_force_ptr(es,system_name);
          }

          if (system.f_epsilon_ptr!=NULL)
          {
              epsilon = system.f_epsilon_ptr(es,system_name);
          }

          velu = velocity(0);
          velv = velocity(1);


         

          // SUPG FORMULATION

          //Calculate supg constant
          NumberVectorValue Usupg = velocity;
          /*Real auxvar =0;
          for (unsigned int l=0; l<phi.size(); l++)
          {
              auxvar += fabs( Usupg*dphi[l][qp]);
          }
          Real hugn;
          if (auxvar ==0)
              hugn = elem->hmin()*2/sqrt(libMesh::pi);
          else
              hugn = 2*Usupg.size()/auxvar;

          //hugn = elem->hmin()*2/sqrt(libMesh::pi);
          const Real kinematic_visc = epsilon;
          Real Reugn = Usupg.size()*hugn/2/kinematic_visc;
          Real zeta = (Reugn<=3) ? Reugn/3 : 1;
          Real tsupg = 1/sqrt(pow(2*Usupg.size()/hugn,2)+pow(4*kinematic_visc/(hugn*hugn),2)+4/dt/dt)*0.35;//std=0.35
          //tsupg is apparently optimized for dt=0.1
          Real tlsic = hugn/2*Usupg.size()*zeta;
          */
          const Real hugn = elem->hmin()*2/sqrt(libMesh::pi);
          Real tsupg = 1/sqrt(pow(2*Usupg.size()/hugn,2)+1*pow(4*epsilon/(hugn*hugn),2)+4/dt/dt)*1;

         


           for (unsigned int i=0;i<phi.size();i++)
          {
            Fe(i) += JxW[qp]*(
                                +uold*phi[i][qp]
                                +dt*(-1 + theta)*(duolddx*velu + duolddy*velv)*phi[i][qp]
                                +dt*(duolddx*(dphi[i][qp](0)) + duolddy*(dphi[i][qp](1)))*epsilon*(-1 + theta)
                                +f*phi[i][qp]*dt
                                +tsupg*((dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)*
             (uold + dt*(f + d2uolddx2*(epsilon - epsilon*theta) +
                d2uolddy2*(epsilon - epsilon*theta) + (-1 + theta)*
                 (duolddx*velu + duolddy*velv)))
                               );

             for (unsigned int j=0; j<phi.size(); j++)
             {
             Ke(i,j) += JxW[qp]*(
                                +phi[j][qp]*phi[i][qp]
                                +dt*theta*((dphi[j][qp](0))*velu + (dphi[j][qp](1))*velv)*phi[i][qp]
                                +dt*((dphi[j][qp](0))*(dphi[i][qp](0)) + (dphi[j][qp](1))*(dphi[i][qp](1)))*epsilon*theta
                                +0
                                +tsupg*((dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)*
             (phi[j][qp] + dt*theta*(-((d2phidx2[j][qp] + d2phidy2[j][qp])*epsilon) + (dphi[j][qp](0))*velu +
                (dphi[j][qp](1))*velv))
                               );

             }

            }
          
          /*
          //getting ready for masud
          const Point& qpoint = qpoints[qp];
          const Real b1 = system.bubble1(qpoint);
          const Real b2 = system.bubble2(qpoint);
          const RealGradient grad_b1 = system.grad_bubble1(qpoint);
          const RealGradient grad_b2 = system.grad_bubble2(qpoint);
          Real tchiu = 1/(b2*(velocity*grad_b1)+epsilon*grad_b2*grad_b1);
          Real tmasud = b1*masud_int*tchiu;
          std::cout << tchiu << "  ";




            for (unsigned int i=0;i<phi.size();i++)
            {
            Fe(i) += JxW[qp]*(
                                +uold*phi[i][qp]
                                +dt*(-1 + theta)*(duolddx*velu + duolddy*velv)*phi[i][qp]
                                +dt*(duolddx*(dphi[i][qp](0)) + duolddy*(dphi[i][qp](1)))*epsilon*(-1 + theta)
                                +dt*f*phi[i][qp]
                                +tmasud*uold*((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon + (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)
                                +dt*(-1 + theta)*tmasud*(duolddx*velu + duolddy*velv)*
             ((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon + (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)
                                +-((d2uolddx2 + d2uolddy2)*dt*epsilon*(-1 + theta)*tmasud*
              ((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon + (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv))
                                +dt*f*tmasud*((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon + (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)
                               );

             for (unsigned int j=0; j<phi.size(); j++)
             {
             Ke(i,j) += JxW[qp]*(
                                +phi[j][qp]*phi[i][qp]
                                +dt*theta*((dphi[j][qp](0))*velu + (dphi[j][qp](1))*velv)*phi[i][qp]
                                +dt*((dphi[j][qp](0))*(dphi[i][qp](0)) + (dphi[j][qp](1))*(dphi[i][qp](1)))*epsilon*theta
                                +0
                                +tmasud*phi[j][qp]*((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon + (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)
                                +dt*theta*tmasud*((dphi[j][qp](0))*velu + (dphi[j][qp](1))*velv)*((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon +
              (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv)
                                +-((d2phidx2[j][qp] + d2phidy2[j][qp])*dt*epsilon*theta*tmasud*
              ((d2phidx2[i][qp] + d2phidy2[i][qp])*epsilon + (dphi[i][qp](0))*velu + (dphi[i][qp](1))*velv))
                                +0
                               );

             }

            }*/






          // Now compute the element matrix and RHS contributions.
          //USING THE STREAMLINE DIFFUSION METHOD + CRANK NICHOLSON DISCRETIZATION SCHEME
          //SET delta_supg_lvlset = 0 for standard galerkin //check Claes Jhonson p 184 for details
          /*for (unsigned int i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*(
                                // Mass matrix term (phi[i][qp] + (tau/mod_vel)*(dphi[i][qp]*velocity))
                                //u_old*(phi[i][qp]+delta_supg_lvlset*(velocity*dphi[i][qp])) +
                                u_old*(phi[i][qp] + (ksupg)*(dphi[i][qp]*velocity)) +
                                -.5*dt*(
                                        // Convection term
                                        // (grad_u_old may be complex, so the
                                        // order here is important!)
                                        (grad_u_old*velocity)*(phi[i][qp]+(ksupg)*(velocity*dphi[i][qp]))+

                                        // Diffusion term                   (-)?
                                        epsilon*(grad_u_old*dphi[i][qp])-epsilon*(ksupg)*(dphi[i][qp]*velocity)*(d2udx2_old+d2udy2_old)

                                        )
                                );

              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(
                                      // Mass-matrix
                                      (phi[i][qp]+(ksupg)*(velocity*dphi[i][qp]))*phi[j][qp] +

                                      0.5*dt*(
                                             // Convection term
                                             (velocity*dphi[j][qp])*phi[i][qp] +(ksupg)*(dphi[j][qp]*velocity)*(dphi[i][qp]*velocity)+

                                             // Diffusion term                (-)?
                                             epsilon*(dphi[i][qp]*dphi[j][qp])-epsilon*(ksupg)*(dphi[i][qp]*velocity)*(d2phidx2[j][qp]+d2phidy2[j][qp])

                                             )
                                      );
                }
            }*/
        }

      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method.
      //
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      {
        // The penalty value.
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        /*for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == NULL)
            {
              fe_face->reinit(elem,s);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  const Number value = 0;

                  // RHS contribution
                  for (unsigned int i=0; i<psi.size(); i++)
                    Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];

                  // Matrix contribution
                  for (unsigned int i=0; i<psi.size(); i++)
                    for (unsigned int j=0; j<psi.size(); j++)
                      Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                }
            }*/
      }

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p SparseMatrix::add_matrix()
      // and \p NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  // That concludes the system matrix assembly routine.
}

template class AdvectionDiffusionSystem<LinearImplicitSystem>;
template class AdvectionDiffusionSystem<TransientLinearImplicitSystem>;
