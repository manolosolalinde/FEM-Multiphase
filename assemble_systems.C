
#include "assemble_systems.h"

RealVectorValue normalx_velocity (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      Real h = system.elem_ptr->hmax();
      RealVectorValue velocity(0,0);
      Gradient grad_sgnlevel = *(system.data_var_map_gradient["grad_SgnLevel"]);
      Real lvlset_reinit = *(system.data_var_map["SgnLevel"]);
      Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");
      Real delta_min = interface_obj_ptr->getDeltax_min();

      velocity = 0.2*(interface_obj_ptr->S(lvlset_reinit,delta_min))*grad_sgnlevel;//1*(interface_obj_ptr->S(lvlset_reinit,2*delta_min))
	

      return velocity;
}

Real normalx_force (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      Gradient grad_sgnlevel = *(system.data_var_map_gradient["grad_SgnLevel"]);

      Real force = grad_sgnlevel(0)/system.data.dt;

      return force;
}

RealVectorValue normaly_velocity (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);


      Real h = system.elem_ptr->hmax();
      RealVectorValue velocity(0,0);
      Gradient grad_sgnlevel = *(system.data_var_map_gradient["grad_SgnLevel"]);
      Real lvlset_reinit = *(system.data_var_map["SgnLevel"]);
      Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");
      Real delta_min = interface_obj_ptr->getDeltax_min();

      velocity = 0.2*(interface_obj_ptr->S(lvlset_reinit,delta_min))*grad_sgnlevel;

      return velocity;
}

Real normaly_force (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      Gradient grad_sgnlevel = *(system.data_var_map_gradient["grad_SgnLevel"]);

      Real force = grad_sgnlevel(1)/system.data.dt;

      return force;
}



RealVectorValue reinit_velocity (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      RealVectorValue velocity;
      Real h = system.elem_ptr->hmax();
      //Real lvlset = *(system.data_var_map["Level"]);
      Real lvlset_reinit = *(system.data_var_map["uold"]);
      // OR Real lvlset = system.get_variable("Level");
      Gradient grad_u_old = *(system.data_var_map_gradient["grad_uold"]);

      Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");
      Real delta_min = interface_obj_ptr->getDeltax_min();


      if (grad_u_old.size()==0)
        velocity = 0.2*(interface_obj_ptr->S(lvlset_reinit,delta_min))*(grad_u_old);//may be better w sgn(u_0)
      else
        velocity = 0.2*(interface_obj_ptr->S(lvlset_reinit,delta_min))*(grad_u_old/grad_u_old.size());//may be better w sgn(u_0)

//      if (grad_u_old.size()<0.6)
//          velocity = velocity*0.95;

      return velocity;
}

Real reinit_force (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      //Real lvlset = *(system.data_var_map["Level"]);
      Real lvlset_reinit = *(system.data_var_map["uold"]);
      Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");
      Real delta_min = interface_obj_ptr->getDeltax_min();
      Real h = system.elem_ptr->hmax();
      Real force = 0.2*(interface_obj_ptr->S(lvlset_reinit,delta_min));
      //Gradient mean_grad_u_old = system.data_var_map_mean_gradient["uold"];

//      if (mean_grad_u_old.size()<0.6)
//          force = force*0.95;

//      if (lvlset < -0.2)
//      {
//          Real k = 0.2*(1-(fabs(lvlset) - 0.2)/0.2);
//          if (k<0)
//              k=0;
//          force = k*(interface_obj_ptr->S(lvlset_reinit,2*delta_min));
//
//      }

      return force;
}

Real reinit_epsilon (EquationSystems& es,
              const std::string& system_name)
{
      // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      
      Number epsilon = es.parameters.get<Real>("Artificial Diffusion reinit");

      Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");
//      Number velocity_low_limit = es.parameters.get<Real>("velocity_low_limit");
//      Number diffusion_factor_lowveloc = es.parameters.get<Real>("diffusion_factor_lowveloc");
//
      Real lvlset = *(system.data_var_map["Level"]);
      Real delta_min = interface_obj_ptr->getDeltax_min();
      //RealVectorValue velocity = system.data_var_map_mean_gradient["uold"];

//      if (velocity.size() <  velocity_low_limit)
//      {
//          //epsilon=epsilon*diffusion_factor_lowveloc;
//      }

      //Elem* el = system.get_element_ptr();
//      if(lvlset< -delta_min)
//          epsilon=epsilon*10;
      

      

      return epsilon;
}


RealVectorValue lvlset_velocity (EquationSystems& es,
              const std::string& system_name)
{
          // Get a reference to the system object.
      TransientAdvectionDiffusionSystem & system =
                es.get_system<TransientAdvectionDiffusionSystem>(system_name);

      Real vel_u = *(system.data_var_map["u"]);
      Real vel_v = *(system.data_var_map["v"]);

    RealVectorValue velocity = RealVectorValue (vel_u,vel_v);
    return velocity;
}


void assemble_stokes (EquationSystems& es,
                      const std::string& system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert (system_name == "Navier-Stokes");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("Navier-Stokes");
  LinearImplicitSystem & system_heaviside =
    es.get_system<LinearImplicitSystem> ("Heaviside Function");
  LinearImplicitSystem & system_kurvature =
    es.get_system<LinearImplicitSystem> ("Kurvature");
  LinearImplicitSystem & system_normalx =
          es.get_system<LinearImplicitSystem> ("normalx SYSTEM");
  LinearImplicitSystem & system_normaly =
          es.get_system<LinearImplicitSystem> ("normaly SYSTEM");
  LinearImplicitSystem & system_reinit =
          es.get_system<LinearImplicitSystem> ("Reinit_LevelSet");
//  LinearImplicitSystem & system_wforce =
//          es.get_system<LinearImplicitSystem> ("Whole Force System");


//  const unsigned int wforcex_var = system_wforce.variable_number ("wforcex");
//  const unsigned int wforcey_var = system_wforce.variable_number ("wforcey");
//  const unsigned int wforcez_var = (dim==3) ? system_wforce.variable_number ("wforcez") : 0;


  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = navier_stokes_system.variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  AutoPtr<FEBase> fe_vel_face  (FEBase::build(dim, fe_vel_type));

  // Build a Finite Element object of the specified type for
  // the pressure variables.
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  QGauss qface (dim-1, fe_vel_type.default_quadrature_order());

  fe_vel_face->attach_quadrature_rule (&qface);
  fe_pres_face->attach_quadrature_rule (&qface);

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe_vel->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<Real> >& d2phidx2 = fe_vel->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe_vel->get_d2phidy2();
  const std::vector<std::vector<Real> >& d2phidxdy = fe_vel->get_d2phidxdy();
  const std::vector<std::vector<Real> >& dphidx = fe_vel->get_dphidx();
  const std::vector<std::vector<Real> >& dphidy = fe_vel->get_dphidy();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();


  const std::vector<Real>& JxW_vel_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();
  const std::vector<Real>& JxW_pres_face = fe_pres_face->get_JxW();
  const std::vector<std::vector<Real> >& psi_face = fe_pres_face->get_phi();
  const std::vector<Point>& qface_points = fe_vel_face->get_xyz();



  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = navier_stokes_system.get_dof_map();
  const DofMap & dof_map_6 = system_heaviside.get_dof_map();
  const DofMap & dof_map_7 = system_kurvature.get_dof_map();
  const DofMap & dof_map_8 = system_normalx.get_dof_map();
  const DofMap & dof_map_9 = system_normaly.get_dof_map();
  const DofMap & dof_map_10 = system_reinit.get_dof_map();
//  const DofMap & dof_map_wforce = system_wforce.get_dof_map();
  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_p;
//  std::vector<unsigned int> dof_indices_wforcex;
//  std::vector<unsigned int> dof_indices_wforcey;
//  std::vector<unsigned int> dof_indices_wforcez;

  std::vector<unsigned int> dof_indices_heaviside;
  std::vector<unsigned int> dof_indices_kurvature;
  std::vector<unsigned int> dof_indices_normalx;
  std::vector<unsigned int> dof_indices_normaly;
  std::vector<unsigned int> dof_indices_reinit;


  const Real dt   = es.parameters.get<Real>("dt");
  Interface * iop = es.parameters.get<Interface*>("interface_obj_ptr");
  Real maxcfl = 0;
  Real maxtsupg = 0;
  // const Real time  = es.parameters.get<Real>("time");
  const Real theta = 1;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      dof_map_6.dof_indices (elem, dof_indices_heaviside);
      dof_map_7.dof_indices (elem, dof_indices_kurvature);
      dof_map_8.dof_indices(elem,dof_indices_normalx);
      dof_map_9.dof_indices(elem,dof_indices_normaly);
      dof_map_10.dof_indices(elem,dof_indices_reinit);

//      dof_map_wforce.dof_indices (elem, dof_indices_wforcex, wforcex_var);
//      dof_map_wforce.dof_indices (elem, dof_indices_wforcey, wforcey_var);
//      if (dim==3)
//          dof_map_wforce.dof_indices (elem, dof_indices_wforcez, wforcez_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      const unsigned int n_kurvature_dofs = dof_indices_kurvature.size();
      const unsigned int n_normalx_dofs = dof_indices_normalx.size();
      const unsigned int n_normaly_dofs = dof_indices_normaly.size();

      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);


      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);


      Number Fr = es.parameters.get<Real>("Froude");
      Number Re = es.parameters.get<Real>("Reynolds");
      Number We = es.parameters.get<Real>("Weber");
      Number gravity_on = es.parameters.get<Number>("Gravity On");
      Number stokes_tau_factor = es.parameters.get<Real>("stokes_tau_factor");


      //SUPG CONSTANT CALCULATION
      /*RealVectorValue integral_veloc;
      integral_veloc.zero();
      for (unsigned int qp=0;qp<qrule.n_points();qp++)
      {
          Real   u =0;
          Real   v =0;
          Real   uold =0;
          Real   vold =0;
          for (unsigned int l=0; l<n_u_dofs ; l++)
          {
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              uold += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              vold += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
          }
          NumberVectorValue qp_velocity (uold,     vold);
          integral_veloc+=JxW[qp]*qp_velocity;
      }

      RealVectorValue supg_velocity = integral_veloc/elem->volume();
      
      Real tau, Pe_x, Pe_y, alfa_x, alfa_y, mod_vel;
      Real del_x;
      Real vel_u = supg_velocity(0);
      Real vel_v = supg_velocity(1);

      del_x = elem->hmin(); //v_f = vel_u
      Pe_x=(vel_u*del_x)/(2*1/1/Re);
      Pe_y=(vel_v*del_x)/(2*1/1/Re);
      alfa_x=(1/tanh(Pe_x))-1/Pe_x;
      alfa_y=(1/tanh(Pe_y))-1/Pe_y;
      if (fabs(Pe_x) < 1e-9)
          alfa_x=0;
      if (fabs(Pe_y) < 1e-9)
          alfa_y=0;
      tau=(alfa_x*vel_u*del_x + alfa_y*vel_v*del_x)/2*stokes_tau_factor;
      mod_vel=supg_velocity.size();
      Real tsupg=tau/(mod_vel*mod_vel);
      if (mod_vel==0)
            tsupg=0;
      else
          tsupg = tau/(mod_vel*mod_vel);



          Real auxvar =0;
          for (unsigned int l=0; l<n_u_dofs; l++)
          {
              auxvar += abs( U*dphi[l][qp]);
          }
          Real hugn = 2*U.size()/auxvar;
          const Real hugn = elem->hmin()*2/sqrt(libMesh::pi);
          const Real kinematic_visc = Mu/Rho/Re;
          Real Reugn = U.size()*hugn/2/kinematic_visc;
          Real zeta = (Reugn<=3) ? Reugn/3 : 1;
          Real tsupg = 1/sqrt(pow(2*U.size()/hugn,2)+1*pow(4*kinematic_visc/(hugn*hugn),2)+4/dt/dt)*stokes_tau_factor;
          Real tlsic = hugn/2*U.size()*zeta;*/


      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the solution & its gradient at the previous timestep.
          Number   u = 0., uold = 0.;
          Number   v = 0., vold = 0.;
          Number   p =0.,pold = 0.;
          Number RhoA =es.parameters.get<Real>("RhoA"); //Density
          Number RhoB =es.parameters.get<Real>("RhoB");; //Density
          Number MuA = es.parameters.get<Real>("MuA");; //Viscosity
          Number MuB = es.parameters.get<Real>("MuB");; //Viscosity
          Number kurvature = 0;
          Gradient grad_kurvature;
          Number normal_x =0;
          Number normal_y =0;
          Number lvlset = 0;
          Gradient grad_H;
          Number heaviside =0;
          Gradient grad_u, grad_u_old;
          Gradient grad_v, grad_v_old;
          Gradient grad_p, grad_p_old;
          Number d2udx2 = 0;
          Number d2udy2 = 0;
          Number d2vdx2 = 0;
          Number d2vdy2 = 0;
          Number d2uolddx2 = 0;
          Number d2uolddy2 = 0;
          Number d2volddx2 = 0;
          Number d2volddy2 = 0;
//          Real wforcex = 0;
//          Real wforcey = 0;
//          Real wforcez = 0;


          // Compute the velocity & its gradient from the previous timestep
          // and the old Newton iterate.
          for (unsigned int l=0; l<n_u_dofs; l++)
            {
              // From the old timestep:

              uold += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              vold += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
              d2uolddx2 += d2phidx2[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              d2uolddy2 += d2phidy2[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              d2volddx2 += d2phidx2[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
              d2volddy2 += d2phidy2[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
              grad_u_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_u[l]));
              grad_v_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_v[l]));

              // From the previous Newton iterate:
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              d2udx2 += d2phidx2[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              d2udy2 += d2phidy2[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              d2vdx2 += d2phidx2[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              d2vdy2 += d2phidy2[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));
            }

          for(unsigned int l=0; l<n_kurvature_dofs;l++)
          {
              kurvature += phi[l][qp]*system_kurvature.current_solution(dof_indices_kurvature[l]);
              grad_kurvature.add_scaled (dphi[l][qp],system_kurvature.current_solution (dof_indices_kurvature[l]));

          }


          for (unsigned int l=0; l<n_normalx_dofs;l++)
          {
              normal_x += phi[l][qp]*system_normalx.current_solution(dof_indices_normalx[l]);
          }
          for (unsigned int l=0; l<n_normaly_dofs;l++)
          {
              normal_y += phi[l][qp]*system_normaly.current_solution(dof_indices_normaly[l]);
          }
          for (unsigned int l=0; l<n_normalx_dofs;l++)
          {
              lvlset += phi[l][qp]*system_reinit.current_solution(dof_indices_reinit[l]);
          }


          grad_H.zero();

          for(unsigned int l=0; l<dof_indices_heaviside.size(); l++)
            {
                heaviside += psi[l][qp]*system_heaviside.current_solution(dof_indices_heaviside[l]);
                grad_H.add_scaled (dpsi[l][qp],system_heaviside.current_solution (dof_indices_heaviside[l]));
            }
          for (unsigned int l=0; l<n_p_dofs; l++)
            {
              pold += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);
              p += psi[l][qp]*navier_stokes_system.current_solution (dof_indices_p[l]);
              grad_p.add_scaled (dpsi[l][qp],navier_stokes_system.current_solution (dof_indices_p[l]));
              grad_p_old.add_scaled (dpsi[l][qp],navier_stokes_system.old_solution (dof_indices_p[l]));
            }
//          for (unsigned int l=0; l<n_p_dofs;l++)
//          {
//              wforcex += psi[l][qp]*system_wforce.current_solution (dof_indices_wforcex[l]);
//              wforcey += psi[l][qp]*system_wforce.current_solution (dof_indices_wforcey[l]);
//              if(dim==3)
//                  wforcez += psi[l][qp]*system_wforce.current_solution (dof_indices_wforcez[l]);
//          }


          // Definitions for convenience.  It is sometimes simpler to do a
          // dot product if you have the full vector at your disposal.
          const NumberVectorValue U_old (uold, vold);
          const NumberVectorValue U     (u,     v);
          const Number  dudx = grad_u(0);
          const Number  dudy = grad_u(1);
          const Number  dvdx = grad_v(0);
          const Number  dvdy = grad_v(1);
          const Number  dpdx = grad_p(0);
          const Number  dpdy = grad_p(1);
          const Number  duolddx = grad_u_old(0);
          const Number  duolddy = grad_u_old(1);
          const Number  dvolddx = grad_v_old(0);
          const Number  dvolddy = grad_v_old(1);
          const Number  dpolddx = grad_p_old(0);
          const Number  dpolddy = grad_p_old(1);

          //Density and viscosity Calculation
          heaviside = iop->H(lvlset); //FIXME heaviside overwrite
          Number Rho = RhoB + (RhoA - RhoB)*heaviside;
          Number Mu = MuB + (MuA - MuB)*heaviside;



          //Calculate supg constant
          NumberVectorValue Usupg = U;
          /*const Real hugn = elem->hmax()*2/sqrt(libMesh::pi);
          const Real kinematic_visc = Mu/Rho/Re;
          Real Reugn = Usupg.size()*hugn/2/kinematic_visc;
          Real zeta = (Reugn<=3) ? Reugn/3 : 1;
          Real tsupg = 1/sqrt(pow(2*Usupg.size()/hugn,2)+pow(4*kinematic_visc/(hugn*hugn),2)+4/dt/dt);
          Real tlsic = hugn/2*Usupg.size()*zeta;*/
          Real auxvar =0;
          for (unsigned int l=0; l<n_u_dofs; l++)
          {
              auxvar += fabs( Usupg*dphi[l][qp]);
          }
          Real hugn;
          if (auxvar ==0)
              hugn = elem->hmin()*2/sqrt(libMesh::pi);
          else
              hugn = 2*Usupg.size()/auxvar;
          const Real kinematic_visc = Mu/Rho/Re;
          Real Reugn = Usupg.size()*hugn/2/kinematic_visc;
          Real zeta = (Reugn<=3) ? Reugn/3 : 1;
          Real tsupg = 1/sqrt(pow(2*Usupg.size()/hugn,2)+pow(4*kinematic_visc/(hugn*hugn),2)+4/dt/dt);
          Real tlsic = hugn/2*Usupg.size()*zeta;

          //Calculate CFL
          Real cfl = fabs(u*dt/elem->hmin())+fabs(v*dt/elem->hmin());
          if (maxcfl<cfl)
              maxcfl = cfl;

          if (maxtsupg<tsupg)
              maxtsupg = tsupg;

          //Real dH = iop->dH(lvlset);
          //kurvature=-kurvature;
          Real d_nx = grad_H(0);
          Real d_ny = grad_H(1);
          Real d_nz = grad_H(2);
          RealVectorValue seminormal(d_nx,d_ny,d_nz);
          RealVectorValue norm_ort(-d_ny,d_nx);
          Real tifce = 0.01;
          Real diffdh = 0.001*grad_H.size();
          //tsupg = 0;


for (unsigned int i=0;i<n_u_dofs;i++)
{
Fu(i) += JxW[qp]*(
                    +(Rho*uold*phi[i][qp])/dt
                    +0
                    +Rho*(dudx*u + dudy*v)*phi[i][qp]
                    +0
                    +(-(1./We)*kurvature*grad_H(0))*phi[i][qp]
                    //+(-(1./We)*wforcex)*phi[i][qp]
                    +(Rho*tsupg*uold*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/dt
                    +0
                    +Rho*tsupg*(dudx*u + dudy*v)*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                    +(-(1./We)*kurvature*grad_H(0))*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                   );

Fv(i) += JxW[qp]*(
                    +(Rho*vold*phi[i][qp])/dt
                    +0
                    +Rho*(dvdx*u + dvdy*v)*phi[i][qp]
                    +0
                    +(-(1./We)*kurvature*grad_H(1) + Rho*(-1./(Fr*Fr))*gravity_on)*phi[i][qp]
                    //+(-(1./We)*wforcey + Rho*(-1./(Fr*Fr))*gravity_on)*phi[i][qp]
                    +(Rho*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)*vold)/dt
                    +0
                    +Rho*tsupg*(dvdx*u + dvdy*v)*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                    +(-(1./We)*kurvature*grad_H(1) + Rho*(-1./(Fr*Fr))*gravity_on)*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                   );

 for (unsigned int j=0; j<n_u_dofs; j++)
 {
 Kuu(i,j) += JxW[qp]*(
                    +(Rho*phi[j][qp]*phi[i][qp])/dt
                    +((2*(dphi[j][qp](0))*(dphi[i][qp](0)) + (dphi[j][qp](1))*(dphi[i][qp](1)))*Mu)/Re
                    +Rho*((dphi[j][qp](0))*u + dudx*phi[j][qp] + (dphi[j][qp](1))*v)*phi[i][qp]
                    +0
                    +0
                    +(Rho*tsupg*phi[j][qp]*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/dt
                    +-(((2*d2phidx2[j][qp] + d2phidy2[j][qp])*Mu*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/Re)
                    +Rho*tsupg*((dphi[j][qp](0))*u + dudx*phi[j][qp] + (dphi[j][qp](1))*v)*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                    +0
                    +tifce*(norm_ort*dphi[j][qp])*(norm_ort*dphi[i][qp])
                    +diffdh*(dphi[j][qp]*dphi[i][qp])
                   );

 Kuv(i,j) += JxW[qp]*(
                    +0
                    +((dphi[j][qp](0))*(dphi[i][qp](1))*Mu)/Re
                    +dudy*Rho*phi[j][qp]*phi[i][qp]
                    +0
                    +0
                    +0
                    +-((d2phidxdy[j][qp]*Mu*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/Re)
                    +dudy*Rho*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)*phi[j][qp]
                    +0
                   );

 Kvu(i,j) += JxW[qp]*(
                    +0
                    +((dphi[j][qp](1))*(dphi[i][qp](0))*Mu)/Re
                    +dvdx*Rho*phi[j][qp]*phi[i][qp]
                    +0
                    +0
                    +0
                    +-((d2phidxdy[j][qp]*Mu*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/Re)
                    +dvdx*Rho*tsupg*phi[j][qp]*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                    +0
                   );

 Kvv(i,j) += JxW[qp]*(
                    +(Rho*phi[j][qp]*phi[i][qp])/dt
                    +(((dphi[j][qp](0))*(dphi[i][qp](0)) + 2*(dphi[j][qp](1))*(dphi[i][qp](1)))*Mu)/Re
                    +Rho*((dphi[j][qp](0))*u + (dphi[j][qp](1))*v + dvdy*phi[j][qp])*phi[i][qp]
                    +0
                    +0
                    +(Rho*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)*phi[j][qp])/dt
                    +-(((d2phidx2[j][qp] + 2*d2phidy2[j][qp])*Mu*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/Re)
                    +Rho*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)*((dphi[j][qp](0))*u + (dphi[j][qp](1))*v + dvdy*phi[j][qp])
                    +0
                    +tifce*(norm_ort*dphi[j][qp])*(norm_ort*dphi[i][qp])
                    +diffdh*(dphi[j][qp]*dphi[i][qp])
                   );

 }

 for (unsigned int j=0;j<n_p_dofs;j++)
 {
 Kup(i,j) += JxW[qp]*(
                    +0
                    +-((dphi[i][qp](0))*psi[j][qp])
                    +0
                    +0
                    +0
                    +0
                    +(dpsi[j][qp](0))*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                    +0
                    +0
                   );

 Kvp(i,j) += JxW[qp]*(
                    +0
                    +-((dphi[i][qp](1))*psi[j][qp])
                    +0
                    +0
                    +0
                    +0
                    +(dpsi[j][qp](1))*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                    +0
                    +0
                   );

 }
}

for (unsigned int i=0;i<n_p_dofs;i++)
{
 for (unsigned int j=0;j<n_u_dofs;j++)
 {
 Kpu(i,j) += JxW[qp]*(
                    +0
                    +0
                    +0
                    +(dphi[j][qp](0))*psi[i][qp]
                    +0
                    +0
                    +0
                    +0
                    +0
                   );

 Kpv(i,j) += JxW[qp]*(
                    +0
                    +0
                    +0
                    +(dphi[j][qp](1))*psi[i][qp]
                    +0
                    +0
                    +0
                    +0
                    +0
                   );

 }
}




        } // end of the quadrature point qp-loop

      //TODO Boundaries genericos
      RealVectorValue (*stokes_dv_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< RealVectorValue (*)(const Real , const Real , const Real)>("stokes_dv_boundary_fptr");
      Real (*stokes_dp_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< Real (*)(const Real , const Real , const Real)>("stokes_dp_boundary_fptr");
      bool (*stokes_is_dv_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< bool (*)(const Real , const Real , const Real)>("stokes_is_dv_boundary_fptr");
      bool (*stokes_is_dp_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< bool (*)(const Real , const Real , const Real)>("stokes_is_dp_boundary_fptr");

      //Begin boundary conditions
      {
          
        const Real penalty = 1.e10;
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == NULL)
            {

              AutoPtr<Elem> side (elem->build_side(s));

              // Loop over the nodes on the side.
              for (unsigned int ns=0; ns< side->n_nodes(); ns++)
                {

                  // Loop over the nodes of the element
                  for (unsigned int n=0; n<elem->n_nodes(); n++)
                  {
                    if (elem->node(n) == side->node(ns))
                    {

                        Node* sidenode = elem->get_node(n);
                        Real x = sidenode->operator ()(0);
                        Real y = sidenode->operator ()(1);
                        Real z = sidenode->operator ()(2);

                       if (stokes_dv_boundary_fptr!=NULL && stokes_is_dv_boundary_fptr!=NULL)
                       {
                        
                            if (n<n_v_dofs && (*stokes_is_dv_boundary_fptr)(x,y,z)==true)
                            {
                                const RealVectorValue vector = (*stokes_dv_boundary_fptr)(x,y,z);
                                // Matrix contribution.
                                Kuu(n,n) += penalty;
                                Kvv(n,n) += penalty;

                                // Right-hand-side contribution.
                                Fu(n) += penalty*vector(0);
                                Fv(n) += penalty*vector(1);

                                //TODO 3D
                            }
                       }
                       else
                       {
                               if (n<n_v_dofs)
                               {
                                Kuu(n,n) += penalty;
                                Kvv(n,n) += penalty;

                                //TODO 3D
                               }

                       }

                       if (stokes_dp_boundary_fptr!=NULL && stokes_is_dp_boundary_fptr!=NULL)
                       {
                            if (n<n_p_dofs && (*stokes_is_dp_boundary_fptr)(x,y,z)==true)
                            {
                                const Real pressure = (*stokes_dp_boundary_fptr)(x,y,z);

                                Kpp(n,n) += penalty;
                                Fp(n)    += penalty*pressure;

                            }
                       }
                       else
                       {
                            const unsigned int pressure_node = 0;
                            if (elem->node(n) == pressure_node)
                            {
                              const Real p_value = 0.0;
                              Kpp(n,n) += penalty;
                              Fp(n)    += penalty*p_value;
                            }

                       }

                    }
                    
                  }
                } // end face node loop
            } // end if (elem->neighbor(side) == NULL)

      }
      // end boundary condition section

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p PetscMatrix::add_matrix()
      // and \p PetscVector::add_vector() members do this for us.
      navier_stokes_system.matrix->add_matrix (Ke, dof_indices);
      navier_stokes_system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

  //store max cfl
  es.parameters.set<Real>("CFL")=maxcfl;
  es.parameters.set<Real>("maxtsupg")=maxtsupg;

  return;
}


void assemble_normal (EquationSystems& es, //TODO functions is deprecated
                      const std::string& system_name)
{
/*
    FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  FEType fe_pres_type = navier_stokes_system.variable_type(p_var);

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));

  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
const std::vector<Real>& JxW = fe_vel->get_JxW();
  */
  libmesh_assert (system_name == "LevelSet Normal");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_normal = es.get_system<LinearImplicitSystem> ("LevelSet Normal");

  TransientLinearImplicitSystem & system_reinit =
    es.get_system<TransientLinearImplicitSystem> ("Reinit_LevelSet");


  // Numeric ids corresponding to each variable in the system
  const unsigned int normal_x_var = system_normal.variable_number ("xnormal");
  const unsigned int normal_y_var = system_normal.variable_number ("wynormal");
  const unsigned int kurv_var = system_normal.variable_number ("kurv");
  const unsigned int auxf_var = system_normal.variable_number ("auxf");

  const unsigned int lvlset_var = system_reinit.variable_number ("SgnLevel");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = system_normal.variable_type(normal_x_var);
  FEType fe_pres_type = system_normal.variable_type(auxf_var);



  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));


  QGauss qrule (dim, fe_vel_type.default_quadrature_order());//fe_vel_type.default_quadrature_order()

  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe_vel->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<Real> >& d2phidx2 = fe_vel->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe_vel->get_d2phidy2();

  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();


  const DofMap & dof_map = system_normal.get_dof_map();
  const DofMap & dof_map_4 = system_reinit.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke),Kuk(Ke), Kuf(Ke),
    Kvu(Ke), Kvv(Ke),Kvk(Ke), Kvf(Ke),
    Kku(Ke), Kkv(Ke),Kkk(Ke), Kkf(Ke),
    Kfu(Ke), Kfv(Ke),Kfk(Ke), Kff(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fk(Fe),
    Ff(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_normal_x;
  std::vector<unsigned int> dof_indices_normal_y;
  std::vector<unsigned int> dof_indices_kurv;
  std::vector<unsigned int> dof_indices_auxf;
  std::vector<unsigned int> dof_indices_lvlset;

  Interface * interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");

  Real time = es.parameters.get<Real>("time");
  Real dt = es.parameters.get<Real>("dt");


  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_normal_x, normal_x_var);
      dof_map.dof_indices (elem, dof_indices_normal_y, normal_y_var);
      dof_map.dof_indices (elem, dof_indices_kurv, kurv_var);
      dof_map.dof_indices (elem, dof_indices_auxf, auxf_var);

      dof_map_4.dof_indices (elem, dof_indices_lvlset, lvlset_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_normal_x.size();
      const unsigned int n_v_dofs = dof_indices_normal_y.size();
      const unsigned int n_kurv_dofs = dof_indices_kurv.size();
      const unsigned int n_auxf_dofs = dof_indices_auxf.size();

      const unsigned int n_lvlset_dofs = dof_indices_lvlset.size();

      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (normal_x_var*n_u_dofs, normal_x_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (normal_x_var*n_u_dofs, normal_y_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kuk.reposition (normal_x_var*n_u_dofs, kurv_var*n_u_dofs, n_u_dofs, n_kurv_dofs);
      Kuf.reposition (normal_x_var*n_u_dofs, auxf_var*n_u_dofs, n_u_dofs, n_auxf_dofs);

      Kvu.reposition (normal_y_var*n_v_dofs, normal_x_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (normal_y_var*n_v_dofs, normal_y_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvk.reposition (normal_y_var*n_v_dofs, kurv_var*n_v_dofs, n_v_dofs, n_kurv_dofs);
      Kvf.reposition (normal_y_var*n_v_dofs, auxf_var*n_v_dofs, n_v_dofs, n_auxf_dofs);

      Kku.reposition (kurv_var*n_v_dofs, normal_x_var*n_kurv_dofs, n_kurv_dofs, n_u_dofs);
      Kkv.reposition (kurv_var*n_v_dofs, normal_y_var*n_kurv_dofs, n_kurv_dofs, n_v_dofs);
      Kkk.reposition (kurv_var*n_v_dofs, kurv_var*n_kurv_dofs, n_kurv_dofs, n_kurv_dofs);
      Kkf.reposition (kurv_var*n_v_dofs, auxf_var*n_kurv_dofs, n_kurv_dofs, n_auxf_dofs);

      Kfu.reposition (auxf_var*n_u_dofs, normal_x_var*n_auxf_dofs, n_auxf_dofs, n_u_dofs);
      Kfv.reposition (auxf_var*n_u_dofs, normal_y_var*n_auxf_dofs, n_auxf_dofs, n_v_dofs);
      Kfk.reposition (auxf_var*n_u_dofs, kurv_var*n_auxf_dofs, n_auxf_dofs, n_kurv_dofs);
      Kff.reposition (auxf_var*n_u_dofs, auxf_var*n_auxf_dofs, n_auxf_dofs, n_auxf_dofs);


      Fu.reposition (normal_x_var*n_u_dofs, n_u_dofs);
      Fv.reposition (normal_y_var*n_u_dofs, n_v_dofs);
      Fk.reposition (kurv_var*n_u_dofs, n_kurv_dofs);
      Ff.reposition (auxf_var*n_u_dofs, n_auxf_dofs);


      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the solution & its gradient at the previous timestep.

          Number lvlset_value = 0;
          Gradient grad_lvlset;
          grad_lvlset.zero();

          for (unsigned int l=0; l<n_lvlset_dofs; l++)
            {
              grad_lvlset.add_scaled (dphi[l][qp],system_reinit.current_solution (dof_indices_lvlset[l]));
              lvlset_value = phi[l][qp]*system_reinit.current_solution(dof_indices_lvlset[l]);
            }
          NumberVectorValue N_lvlset = grad_lvlset;// or grad_lvlset.unit();
          NumberVectorValue N_lvlset_unit = N_lvlset;
          if (N_lvlset.size()>0)
              N_lvlset_unit.unit();


          const Number normal_x = N_lvlset(0);
          const Number normal_y = N_lvlset(1);
          Real epsilon = interface_obj_ptr->getEpsilon()*5;
          Real dH = interface_obj_ptr->dH(lvlset_value,epsilon);
          dH = dH*epsilon;
          if (dH>1)
              libmesh_error();
          const Real diff1 = 0;//pow(0.3,time/dt+8);//0.000005;//0.0000000005;
          const Real diff2 = 0;
          

          for (unsigned int i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fu(i) += JxW[qp]*(normal_x*phi[i][qp]);
              Fv(i) += JxW[qp]*(normal_y*phi[i][qp]);
              //Fu(i) += -JxW[qp]*(lvlset_value*dphi[i][qp](0));
              //Fv(i) += -JxW[qp]*(lvlset_value*dphi[i][qp](1));

              Fk(i) += JxW[qp]*(0);
              


              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                       -dphi[i][qp]*dphi[j][qp]*diff1*dH
                                       //-(N_lvlset_unit*dphi[i][qp])*(N_lvlset_unit*dphi[j][qp])*diff1
                                        );
                  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                       -dphi[i][qp]*dphi[j][qp]*diff1*dH
                                       //-(N_lvlset_unit*dphi[i][qp])*(N_lvlset_unit*dphi[j][qp])*diff1
                                        );

                  Kkk(i,j) += JxW[qp]*(-dphi[i][qp]*dphi[j][qp]*dH*diff2
                                       + phi[i][qp]*phi[j][qp]);

                  Kku(i,j) += JxW[qp]*(-phi[i][qp]*dphi[j][qp](0));
                  Kkv(i,j) += JxW[qp]*(-phi[i][qp]*dphi[j][qp](1));
                }
              for (unsigned int j=0; j<psi.size(); j++)
                {
                  //Kuf(i,j) += JxW[qp]*(dphi[i][qp](0)*psi[j][qp]);
                  //Kvf(i,j) += JxW[qp]*(dphi[i][qp](1)*psi[j][qp]);
                }
            }

          for (unsigned int i=0; i<phi.size(); i++)
            {

              Ff(i) += JxW[qp]*dH*phi[i][qp];
              
              for (unsigned int j=0; j<phi.size(); j++)
              {
                  //Kfu(i,j) += JxW[qp]*(-psi[i][qp]*dphi[j][qp](0));
                  //Kfv(i,j) += JxW[qp]*(-psi[i][qp]*dphi[j][qp](1));
                  //Kfk(i,j) += JxW[qp]*(psi[i][qp]*phi[j][qp]);
                  //Kfk(i,j) +=JxW[qp]*(psi[i][qp]*phi[j][qp]*dH*0.001);
              }

              for (unsigned int j=0; j<phi.size(); j++)
              {
                  Kff(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
              }
            }




      }

      //boundary conditions
      {

//        const Real penalty = 1.e10;
//
//        for (unsigned int i=0; i<psi.size(); i++)
//        {
//            Real lvlset_value = system_reinit.current_solution(dof_indices_lvlset[i]);
//            Real dH = interface_obj_ptr->dH(lvlset_value);
//            if (dH == 0)
//            {
//                Kff(i,i)+=penalty;
//                Ff(i)   +=penalty*0;
//            }
//        }

//        const bool pin_pressure = true;
//        if (pin_pressure)
//          {
//            const unsigned int pressure_node = 0;
//            const Real f_value               = 0.0;
//            for (unsigned int c=0; c<elem->n_nodes(); c++)
//              if (elem->node(c) == pressure_node)
//                {
//                  Kff(c,c) += penalty;
//                  Ff(c)    += penalty*f_value;
//                }
//          }
      }


      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system_normal.matrix->add_matrix (Ke, dof_indices);
      system_normal.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

  // That's it.
  return;
}


void assemble_heaviside (EquationSystems& es,
                      const std::string& system_name)
{
  libmesh_assert (system_name == "Heaviside Function");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_heaviside = es.get_system<LinearImplicitSystem> ("Heaviside Function");

  TransientLinearImplicitSystem & system_reinit =
    es.get_system<TransientLinearImplicitSystem> ("Reinit_LevelSet");


  // Numeric ids corresponding to each variable in the system
  const unsigned int heaviside_var = system_heaviside.variable_number ("H");
  const unsigned int sgnlevel_var = system_reinit.variable_number("SgnLevel");
  
  // Get the Finite Element type 
  FEType fe_pres_type = system_heaviside.variable_type(heaviside_var);
  FEType fe_vel_type = system_reinit.variable_type(sgnlevel_var);


  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  AutoPtr<FEBase> fe_pres  (FEBase::build(dim, fe_pres_type));


  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe_pres->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<Real> >& d2phidx2 = fe_vel->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe_vel->get_d2phidy2();

  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  const DofMap & dof_map = system_heaviside.get_dof_map();
  const DofMap & dof_map_4 = system_reinit.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_reinit;

  Interface * interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map_4.dof_indices (elem, dof_indices_reinit);

      fe_pres->reinit (elem);
      fe_vel->reinit(elem);

      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
          // Values to hold the solution & its gradient at the previous timestep.

          Number lvlset_value = 0;
          for (unsigned int l=0; l<phi.size(); l++)
            {
              lvlset_value += phi[l][qp]*system_reinit.current_solution  (dof_indices_reinit[l]);
            }

          Real H = interface_obj_ptr->H(lvlset_value);

        for (unsigned int i=0; i<psi.size(); i++)
        {
              Fe(i) += JxW[qp]*psi[i][qp]*H;
            
        
              for (unsigned int j=0; j<psi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*psi[i][qp]*psi[j][qp];
              }
         }
      }


        /*for (unsigned int i=0; i<psi.size(); i++)
        {
	       
            Real lvlset_value = system_reinit.current_solution  (dof_indices_reinit[i]);
            Real H = interface_obj_ptr->H(lvlset_value);

              Fe(i) = H;
              Ke(i,i) = 1;

        }*/
	


      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p PetscMatrix::add_matrix()
      // and \p PetscVector::add_vector() members do this for us.
      system_heaviside.matrix->add_matrix (Ke, dof_indices);
      system_heaviside.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

  // That's it.
  return;
}

void assemble_kurvature (EquationSystems& es,
                  const std::string& system_name)
{
  libmesh_assert (system_name == "Kurvature");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_kurvature =
    es.get_system<LinearImplicitSystem> ("Kurvature");
  TransientAdvectionDiffusionSystem & system_normalx =
    es.get_system<TransientAdvectionDiffusionSystem> ("normalx SYSTEM");
  TransientAdvectionDiffusionSystem & system_normaly =
    es.get_system<TransientAdvectionDiffusionSystem> ("normaly SYSTEM");


  FEType fe_type = system_kurvature.variable_type(0);

  AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));

  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule      (&qrule);

  const std::vector<Real>& JxW      = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<Real> >& dphidx = fe->get_dphidx();
  const std::vector<std::vector<Real> >& dphidy = fe->get_dphidy();
  const std::vector<std::vector<Real> >& d2phidx2 = fe->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe->get_d2phidy2();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<Point>& qp_xyz = fe->get_xyz();

  
  const DofMap& dof_map = system_kurvature.get_dof_map();
  const DofMap & dof_map_nx = system_normalx.get_dof_map();
  const DofMap & dof_map_ny = system_normaly.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_nx;
  std::vector<unsigned int> dof_indices_ny;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map_nx.dof_indices(elem,dof_indices_nx);
      dof_map_ny.dof_indices(elem,dof_indices_ny);

      fe->reinit (elem);

      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          Real nx = 0;
          Real ny = 0;
          Gradient grad_nx;
          Gradient grad_ny;


          for (unsigned int l=0; l<phi.size(); l++)
            {
              nx += phi[l][qp]*system_normalx.current_solution (dof_indices_nx[l]);
              ny += phi[l][qp]*system_normaly.current_solution (dof_indices_ny[l]);
              grad_nx.add_scaled (dphi[l][qp],system_normalx.current_solution (dof_indices_nx[l]));
              grad_ny.add_scaled (dphi[l][qp],system_normaly.current_solution (dof_indices_ny[l]));

            }

          NumberVectorValue ortho_normal = NumberVectorValue(-ny,nx);
          Real tifce = 0.0005;//std=0.001

          for (unsigned int i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*phi[i][qp]*(-ny*(ortho_normal*grad_nx) + nx*(ortho_normal*grad_ny));
              //Fe(i) += JxW[qp]*phi[i][qp]*(grad_nx(0)+grad_ny(1));


              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            //+0.00005*dphi[i][qp]*dphi[j][qp]
                                            +tifce*(ortho_normal*dphi[j][qp])*(ortho_normal*dphi[i][qp])
                                            );
                }
            }
        }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system_kurvature.matrix->add_matrix (Ke, dof_indices);
      system_kurvature.rhs->add_vector    (Fe, dof_indices);
    }

  // That concludes the system matrix assembly routine.
}


void assemble_refinement (EquationSystems& es,
                  const std::string& system_name)
{
  libmesh_assert (system_name == "Refinement");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();


  TransientLinearImplicitSystem & system =
          es.get_system<TransientLinearImplicitSystem>("Refinement");
  TransientLinearImplicitSystem & system_reinit =
    es.get_system<TransientLinearImplicitSystem>("Reinit_LevelSet");

  const DofMap& dof_map = system.get_dof_map();
  const DofMap& dof_map_reinit = system_reinit.get_dof_map();

  FEType fe_vel_type = system.variable_type(0);

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  fe_vel->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();


  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_reinit;

  Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");
  Real local_tol = es.parameters.get<Real>("local tolerance");
  unsigned int max_h_level = es.parameters.get<unsigned int>("max_h_level");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map_reinit.dof_indices (elem, dof_indices_reinit);

      fe_vel->reinit (elem);

      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      Real epsilon = interface_obj_ptr->getEpsilon();

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the old solution & its gradient.
          Number lvlset=0;
	   Gradient grad_heaviside;

          // Compute the old solution & its gradient.
          for (unsigned int l=0; l<phi.size(); l++)
            {
              lvlset   += phi[l][qp]*system_reinit.current_solution  (dof_indices_reinit[l]);
            }

	   Real dH = interface_obj_ptr->dH(lvlset);



          // Now compute the element matrix and RHS contributions.
          for (unsigned int i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*(
                                dH*phi[i][qp]//std dH_ref
                                );


              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(
                                      phi[i][qp]*phi[j][qp]
                                      );
                }
            }

        }//end quadrature point loop
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);


  }//end elem loop

  return;

}

void assemble_corrector (EquationSystems& es, const std::string& system_name)
{


  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system =
          es.get_system<LinearImplicitSystem>("Reinit Corrector");
  TransientLinearImplicitSystem & system_reinit =
    es.get_system<TransientLinearImplicitSystem> ("Reinit_LevelSet");
  TransientLinearImplicitSystem & system_lvlset =
    es.get_system<TransientLinearImplicitSystem> ("LevelSet");

  FEType fe_type = system.variable_type(0);

//  AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
//  QuadratureType quad_type=INVALID_Q_RULE;
//  quad_type = static_cast<QuadratureType>(0);
  // 0 = Gauss-Legendre
//  AutoPtr<QBase> qrule(QBase::build(quad_type, dim, TENTH));
  //QGauss qrule (dim,   fe_type.default_quadrature_order());
//  fe->attach_quadrature_rule      (qrule.get());

  AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule      (&qrule);

  const std::vector<Real>& JxW      = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  const DofMap& dof_map = system.get_dof_map();
  const DofMap& dof_map_lvlset = system_lvlset.get_dof_map();
  const DofMap& dof_map_reinit = system_reinit.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_lvlset;
  std::vector<unsigned int> dof_indices_reinit;

  const Real dt = es.parameters.get<Real>   ("dt_tau");
  const Real time = es.parameters.get<Real> ("tau_time");

  Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map_lvlset.dof_indices (elem, dof_indices_lvlset);
      dof_map_reinit.dof_indices (elem, dof_indices_reinit);

      fe->reinit (elem);

      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      Real numerator =0;
      Real denominator = 0;

      /*for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {

          Number   u_0 =0;
          Number   u =0;
          Number u_old =0;
          Gradient grad_u;

          //OBS: system is system_reinit
          for (unsigned int l=0; l<phi.size(); l++)
          {
              u_0 += phi[l][qp]*system_lvlset.current_solution  (dof_indices_lvlset[l]);
              u      += phi[l][qp]*system_reinit.current_solution  (dof_indices_reinit[l]);
              u_old  += phi[l][qp]*system_reinit.old_solution  (dof_indices_reinit[l]);
              //grad_u_old.add_scaled (dphi[l][qp],system_reinit.old_solution (dof_indices[l]));
              grad_u.add_scaled (dphi[l][qp],system_reinit.current_solution (dof_indices_reinit[l]));
          }
          Real mod_grad = grad_u.size();
          Number dH = interface_obj_ptr->dH(u_old);

          //not sure what is this for
          //if (mod_grad_old>1)
          //    mod_grad_old = 1;

          //numerator += JxW[qp] * dH*(-S0*(mod_grad_lvl-1));
          mod_grad=1; //coment to unchange
          numerator += JxW[qp] * dH* (u - u_old)/dt;
          denominator += JxW[qp] * dH*dH*mod_grad;
      } 


      Real lambda_elem=0;
      if (fabs(denominator) > 0.0000000000000001)
      {
          lambda_elem += (- numerator/denominator);
      }
       else
       {
          lambda_elem = 0;
          if (denominator!=0)
              std::cout << "Lamda_elem has been aproximated to 0 - this may be an error";

       }*/

      // NOTE: The correction that uses lambda element by element introduces noise on the kurvature
          

      
       //Now we'll find the l2 projection of lambda
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {

      Number   u_0 =0;
      Number   u =0;
      Number u_old =0;
      Gradient grad_u;

      for (unsigned int l=0; l<phi.size(); l++)
          {
              u_0 += phi[l][qp]*system_lvlset.current_solution  (dof_indices_lvlset[l]);
              u      += phi[l][qp]*system_reinit.current_solution  (dof_indices_reinit[l]); 
              u_old  += phi[l][qp]*system_reinit.old_solution  (dof_indices_reinit[l]);//
             grad_u.add_scaled (dphi[l][qp],system_reinit.current_solution (dof_indices_reinit[l]));
          }
      Real mod_grad = grad_u.size();

      //not sure what is this for
      //if (mod_grad_old>1)
      //    mod_grad_old = 1;
      mod_grad=1;//coment to unchange
      Real eps = interface_obj_ptr->getEpsilon();
      

      for (unsigned int i=0; i<phi.size(); i++) 
            {
              // The RHS contribution
              //Fe(i) += JxW[qp]*((lambda_elem*interface_obj_ptr->dH(u_old)*dt*mod_grad)*phi[i][qp]
		  Fe(i) += JxW[qp]*(-(u-u_old)*interface_obj_ptr->dH(u_old,eps*1)*(eps*1)*phi[i][qp]
                                );


              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(
                                      phi[i][qp]*phi[j][qp]
				          //+0.00005*dphi[i][qp]*dphi[j][qp]
                                      );
                }
            }

      }//end qp loop

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);


  }//end elem loop

  return;

}

void assemble_test(EquationSystems& es, const std::string& system_name)
{
  libmesh_assert (system_name == "TEST SYSTEM");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_test =
    es.get_system<LinearImplicitSystem> ("TEST SYSTEM");
  TransientAdvectionDiffusionSystem & system_reinit =
    es.get_system<TransientAdvectionDiffusionSystem> ("Reinit_LevelSet");
  TransientAdvectionDiffusionSystem & system_normalx =
    es.get_system<TransientAdvectionDiffusionSystem> ("normalx SYSTEM");
  TransientAdvectionDiffusionSystem & system_normaly =
    es.get_system<TransientAdvectionDiffusionSystem> ("normaly SYSTEM");


  FEType fe_type = system_test.variable_type(0);

  AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));

  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule      (&qrule);

  const std::vector<Real>& JxW      = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<Real> >& dphidx = fe->get_dphidx();
  const std::vector<std::vector<Real> >& dphidy = fe->get_dphidy();
  const std::vector<std::vector<Real> >& d2phidx2 = fe->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe->get_d2phidy2();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<Point>& qp_xyz = fe->get_xyz();


  const DofMap& dof_map = system_test.get_dof_map();
  const DofMap & dof_map_nx = system_normalx.get_dof_map();
  const DofMap & dof_map_ny = system_normaly.get_dof_map();
  const DofMap & dof_map_reinit = system_reinit.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_nx;
  std::vector<unsigned int> dof_indices_ny;
  std::vector<unsigned int> dof_indices_reinit;

  Interface* interface_obj_ptr = es.parameters.get<Interface*>("interface_obj_ptr");


  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map_nx.dof_indices(elem,dof_indices_nx);
      dof_map_ny.dof_indices(elem,dof_indices_ny);
      dof_map_reinit.dof_indices(elem,dof_indices_reinit);

      fe->reinit (elem);

      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          Real nx = 0;
          Real ny = 0;
          Gradient grad_nx;
          Gradient grad_ny;
          Real lvlset = 0;


          for (unsigned int l=0; l<phi.size(); l++)
            {
              nx += phi[l][qp]*system_normalx.current_solution (dof_indices_nx[l]);
              ny += phi[l][qp]*system_normaly.current_solution (dof_indices_ny[l]);
              grad_nx.add_scaled (dphi[l][qp],system_normalx.current_solution (dof_indices_nx[l]));
              grad_ny.add_scaled (dphi[l][qp],system_normaly.current_solution (dof_indices_ny[l]));

            }
          for (unsigned int l=0; l<phi.size(); l++)
            {
              lvlset +=phi[l][qp]*system_reinit.current_solution(dof_indices_reinit[l]);

            }

          Real dH = interface_obj_ptr->dH(lvlset);


          NumberVectorValue ortho_normal = NumberVectorValue(-ny,nx);

          for (unsigned int i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*phi[i][qp]*(dH*nx);


              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
                }
            }
        }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system_test.matrix->add_matrix (Ke, dof_indices);
      system_test.rhs->add_vector    (Fe, dof_indices);
    }

  // That concludes the system matrix assembly routine.
}

void compute_jacobian (const NumericVector<Number>& soln,
                       SparseMatrix<Number>&  jacobian,
                       NonlinearImplicitSystem& sys)
{

    PetscVector<Real> fake_residual;
    compute_navier_stokes(soln,jacobian,fake_residual,sys,2);
  // That's it.
}



void compute_residual (const NumericVector<Number>& soln,
                       NumericVector<Number>& residual,
                       NonlinearImplicitSystem& sys)
{
 // Here we compute the residual R(x) = K(x)*x - f. The current solution
 // x is passed in the soln vector

    //flag = 0 ->implies residual
    PetscMatrix<Real> fake_jacobian;
    compute_navier_stokes(soln,fake_jacobian,residual,sys,1);


}


void compute_navier_stokes (const NumericVector<Number>& soln,
                       SparseMatrix<Number>&  jacobian,
                       NumericVector<Number>& residual,
                       NonlinearImplicitSystem& sys,
                       unsigned int flag)
{

  EquationSystems &es = sys.get_equation_systems();

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  libmesh_assert (dim == 2 || dim==3);
  libmesh_assert (flag==1 || flag==2);

  TransientNonlinearImplicitSystem& navier_stokes_system =
    es.get_system<TransientNonlinearImplicitSystem>("Navier-Stokes");
  LinearImplicitSystem & system_heaviside =
    es.get_system<LinearImplicitSystem> ("Heaviside Function");
  LinearImplicitSystem & system_kurvature =
    es.get_system<LinearImplicitSystem> ("Kurvature");
  TransientLinearImplicitSystem & system_reinit =
    es.get_system<TransientLinearImplicitSystem> ("Reinit_LevelSet");

  LinearImplicitSystem & system_force =
          es.get_system<LinearImplicitSystem> ("Force System");
  LinearImplicitSystem & system_wforce =
          es.get_system<LinearImplicitSystem> ("Whole Force System");
  LinearImplicitSystem & system_jump =
          es.get_system<LinearImplicitSystem> ("Jump Pres System");

  const unsigned int forcex_var = system_force.variable_number ("forcex");
  const unsigned int forcey_var = system_force.variable_number ("forcey");
  //const unsigned int forcez_var = system_force.variable_number ("forcez");
  const unsigned int forcez_var = (dim==3) ? system_force.variable_number ("forcez") : 0;

  const unsigned int wforcex_var = system_wforce.variable_number ("wforcex");
  const unsigned int wforcey_var = system_wforce.variable_number ("wforcey");
  const unsigned int wforcez_var = (dim==3) ? system_wforce.variable_number ("wforcez") : 0;


  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int w_var = (dim==3)?(navier_stokes_system.variable_number ("w")):0;
  const unsigned int p_var = navier_stokes_system.variable_number ("p");

  //const unsigned int force_var = system_force.variable_number ("force");
  //const unsigned int forcey_var = system_force.variable_number ("forcey");
  //const unsigned int forcez_var = (dim==3) ? system_force.variable_number ("forcez") : 0;

  const DofMap& dof_map = navier_stokes_system.get_dof_map();
  const DofMap& dof_map_heaviside = system_heaviside.get_dof_map();
  const DofMap& dof_map_kurvature = system_kurvature.get_dof_map();
  const DofMap& dof_map_reinit  = system_reinit.get_dof_map();
  const DofMap& dof_map_jump = system_jump.get_dof_map();
  const DofMap & dof_map_force = system_force.get_dof_map();
  const DofMap & dof_map_wforce = system_wforce.get_dof_map();


  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);
  FEType fe_pres_type = navier_stokes_system.variable_type(p_var);

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  AutoPtr<FEBase> fe_vel_face  (FEBase::build(dim, fe_vel_type));
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));

  QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  QGauss qface (dim-1, fe_vel_type.default_quadrature_order());

  fe_vel_face->attach_quadrature_rule (&qface);
  fe_pres_face->attach_quadrature_rule (&qface);
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<Real> >& d2phidx2 = fe_vel->get_d2phidx2();
  const std::vector<std::vector<Real> >& d2phidy2 = fe_vel->get_d2phidy2();
  const std::vector<std::vector<Real> >& d2phidxdy = fe_vel->get_d2phidxdy();
  const std::vector<std::vector<Real> >& dphidx = fe_vel->get_dphidx();
  const std::vector<std::vector<Real> >& dphidy = fe_vel->get_dphidy();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();


  const std::vector<Real>& JxW_vel_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();
  const std::vector<Real>& JxW_pres_face = fe_pres_face->get_JxW();
  const std::vector<std::vector<Real> >& psi_face = fe_pres_face->get_phi();
  const std::vector<Point>& qface_points = fe_vel_face->get_xyz();


  DenseVector<Number> Residual;

    DenseSubVector<Number>
    Reu(Residual),
    Rev(Residual),
    Rew(Residual),
    Rep(Residual);


  DenseMatrix<Number> Ke;
    DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke),Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke),Kvp(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke),Kwp(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke),Kpp(Ke);


  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;

  std::vector<unsigned int> dof_indices_jump;
  std::vector<unsigned int> dof_indices_forcex;
  std::vector<unsigned int> dof_indices_forcey;
  std::vector<unsigned int> dof_indices_forcez;
  std::vector<unsigned int> dof_indices_wforcex;
  std::vector<unsigned int> dof_indices_wforcey;
  std::vector<unsigned int> dof_indices_wforcez;

  std::vector<unsigned int> dof_indices_kurvature;
  std::vector<unsigned int> dof_indices_heaviside;
  std::vector<unsigned int> dof_indices_reinit;

  // We will compute the element residual.
  if(flag==1)
      residual.zero();

  const Real dt   = es.parameters.get<Real>("dt");
  Interface * iop = es.parameters.get<Interface*>("interface_obj_ptr");
  Real maxcfl = 0;
  Real maxtsupg = 0;
  // const Real time  = es.parameters.get<Real>("time");
  const Real theta = 1;

  Number RhoA =es.parameters.get<Real>("RhoA"); //Density
  Number RhoB =es.parameters.get<Real>("RhoB");; //Density
  Number MuA = es.parameters.get<Real>("MuA");; //Viscosity
  Number MuB = es.parameters.get<Real>("MuB");; //Viscosity
  Real Fr = es.parameters.get<Real>("Froude");
  Real Re = es.parameters.get<Real>("Reynolds");
  Real We = es.parameters.get<Real>("Weber");
  Real gravity_on = es.parameters.get<Real>("Gravity On");


  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      if(dim==3)
          dof_map.dof_indices (elem, dof_indices_w, w_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      dof_map_heaviside.dof_indices(elem,dof_indices_heaviside);
      dof_map_kurvature.dof_indices(elem,dof_indices_kurvature);
      dof_map_reinit.dof_indices(elem,dof_indices_reinit);

      dof_map_jump.dof_indices(elem,dof_indices_jump);
      dof_map_force.dof_indices (elem, dof_indices_forcex, forcex_var);
      dof_map_force.dof_indices (elem, dof_indices_forcey, forcey_var);
      if (dim==3)
          dof_map_force.dof_indices (elem, dof_indices_forcez, forcez_var);

      dof_map_wforce.dof_indices (elem, dof_indices_wforcex, wforcex_var);
      dof_map_wforce.dof_indices (elem, dof_indices_wforcey, wforcey_var);
      if (dim==3)
          dof_map_wforce.dof_indices (elem, dof_indices_wforcez, wforcez_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      const unsigned int n_kurvature_dofs = dof_indices_kurvature.size();
      const unsigned int n_heaviside_dofs = dof_indices_heaviside.size();
      const unsigned int n_reinit_dofs = dof_indices_reinit.size();


      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      if(flag==2)
      {
          Ke.resize (n_dofs, n_dofs);

          Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
          if(dim==3)
              Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

          Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
          if(dim==3)
              Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kvp.reposition (v_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

          if(dim==3)
          {
          Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kwp.reposition (w_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
          }

          Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
          Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_u_dofs);
          if(dim==3)
              Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_u_dofs);
          Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

      }



      if (flag==1)
      {
          Residual.resize (n_dofs);
          Reu.reposition (u_var*n_u_dofs, n_u_dofs);
          Rev.reposition (v_var*n_u_dofs, n_u_dofs);
          if (dim==3)
              Rew.reposition (w_var*n_u_dofs, n_u_dofs);
          Rep.reposition (p_var*n_u_dofs, n_p_dofs);
      }

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the solution & its gradient at the previous timestep.
          Number   u = 0., uold = 0.;
          Number   v = 0., vold = 0.;
          Number   p =0.,pold = 0.;
          Number kurvature = 0;
          Gradient grad_H;
          Number heaviside =0;
          Gradient grad_u, grad_u_old;
          Gradient grad_v, grad_v_old;
          Gradient grad_p, grad_p_old;
          Number d2udx2 = 0;
          Number d2udy2 = 0;
          Number d2vdx2 = 0;
          Number d2vdy2 = 0;
          Number d2uolddx2 = 0;
          Number d2uolddy2 = 0;
          Number d2volddx2 = 0;
          Number d2volddy2 = 0;
          Number lvlset = 0;
          Gradient normal;
          Number dH = 0;

          Real jump = 0;
          Gradient grad_jump;
          Real forcex =0;
          Real forcey = 0;
          Real forcez = 0;
          Real wforcex =0;
          Real wforcey = 0;
          Real wforcez = 0;

          // Compute the velocity & its gradient from the previous timestep
          // and the old Newton iterate.
          for (unsigned int l=0; l<n_u_dofs; l++)
            {
              // From the old timestep:
              if(flag==1)//just for the residual
              {
                  uold += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
                  vold += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
                  d2uolddx2 += d2phidx2[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
                  d2uolddy2 += d2phidy2[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
                  d2volddx2 += d2phidx2[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
                  d2volddy2 += d2phidy2[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
                  grad_u_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_u[l]));
                  grad_v_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_v[l]));
              }

              // From the previous nonlinear iterate:
              u += phi[l][qp]*soln(dof_indices_u[l]);
              v += phi[l][qp]*soln(dof_indices_v[l]);
              d2udx2 += d2phidx2[l][qp]*soln(dof_indices_u[l]);
              d2udy2 += d2phidy2[l][qp]*soln(dof_indices_u[l]);
              d2vdx2 += d2phidx2[l][qp]*soln(dof_indices_v[l]);
              d2vdy2 += d2phidy2[l][qp]*soln(dof_indices_v[l]);
              grad_u.add_scaled (dphi[l][qp],soln(dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp],soln(dof_indices_v[l]));

            }

          for (unsigned int l=0; l<n_p_dofs; l++)
            {
              pold += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);
              grad_p_old.add_scaled (dpsi[l][qp],navier_stokes_system.old_solution (dof_indices_p[l]));

              // From the previous nonlinear iterate:
              p += psi[l][qp]*soln(dof_indices_p[l]);
              grad_p.add_scaled (dpsi[l][qp],soln(dof_indices_p[l]));
            }

          for (unsigned int l=0; l<n_p_dofs; l++)
          {
              jump += psi[l][qp]*system_jump.current_solution (dof_indices_jump[l]);
              grad_jump.add_scaled (dpsi[l][qp],system_jump.current_solution(dof_indices_jump[l]));
          }


          for (unsigned int l=0; l<n_u_dofs;l++)
          {
              lvlset += phi[l][qp]*system_reinit.current_solution(dof_indices_reinit[l]);
              normal.add_scaled (dphi[l][qp],system_reinit.current_solution (dof_indices_reinit[l]));
              forcex += phi[l][qp]*system_force.current_solution (dof_indices_forcex[l]);
              forcey += phi[l][qp]*system_force.current_solution (dof_indices_forcey[l]);
              if(dim==3)
                  forcez += phi[l][qp]*system_force.current_solution (dof_indices_forcez[l]);
          }

          for (unsigned int l=0; l<n_p_dofs;l++)
          {
              wforcex += psi[l][qp]*system_wforce.current_solution (dof_indices_wforcex[l]);
              wforcey += psi[l][qp]*system_wforce.current_solution (dof_indices_wforcey[l]);
              if(dim==3)
                  wforcez += psi[l][qp]*system_wforce.current_solution (dof_indices_wforcez[l]);
          }

          if (normal.size()!=0)
          {
              normal.unit();
          }
          dH = iop->dH(lvlset);
          Real H = iop->H(lvlset);


          grad_H.zero();
          for(unsigned int l=0; l<n_heaviside_dofs; l++)
            {
                heaviside += psi[l][qp]*system_heaviside.current_solution(dof_indices_heaviside[l]);
                grad_H.add_scaled (dpsi[l][qp],system_heaviside.current_solution (dof_indices_heaviside[l]));
            }

//         if (grad_H.size() > 0.000002467)
              for(unsigned int l=0; l<n_kurvature_dofs;l++)
              {
                  kurvature += phi[l][qp]*system_kurvature.current_solution(dof_indices_kurvature[l]);
              }
          //kurvature = -kurvature;//1/radio;

          // Definitions for convenience.  It is sometimes simpler to do a
          // dot product if you have the full vector at your disposal.
          const NumberVectorValue U_old (uold, vold);
          const NumberVectorValue U     (u,     v);
          const Number  dudx = grad_u(0);
          const Number  dudy = grad_u(1);
          const Number  dvdx = grad_v(0);
          const Number  dvdy = grad_v(1);
          const Number  dpdx = grad_p(0);
          const Number  dpdy = grad_p(1);
          const Number  duolddx = grad_u_old(0);
          const Number  duolddy = grad_u_old(1);
          const Number  dvolddx = grad_v_old(0);
          const Number  dvolddy = grad_v_old(1);
          const Number  dpolddx = grad_p_old(0);
          const Number  dpolddy = grad_p_old(1);

          //Density and viscosity Calculation
          //Number Rho = RhoB + (RhoA - RhoB)*heaviside;
          //Number Mu = MuB + (MuA - MuB)*heaviside;
          Number Rho = RhoB + (RhoA - RhoB)*heaviside;
          Number Mu = MuB + (MuA - MuB)*heaviside;

          //Calculate supg constant
          NumberVectorValue Usupg = U;
          Real auxvar =0;
          for (unsigned int l=0; l<n_u_dofs; l++)
          {
              auxvar += fabs( Usupg*dphi[l][qp]);
          }
          Real hugn;
          if (auxvar ==0)
              hugn = elem->hmin()*2/sqrt(libMesh::pi);
          else
              hugn = 2*Usupg.size()/auxvar;
          const Real kinematic_visc = Mu/Rho/Re;
          const Real dtsupg = dt;
          Real Reugn = Usupg.size()*hugn/2/kinematic_visc;
          Real zeta = (Reugn<=3) ? Reugn/3 : 1;
          Real tsupg = 1/sqrt(pow(2*Usupg.size()/hugn,2)+pow(4*kinematic_visc/(hugn*hugn),2)+4/dtsupg/dtsupg);
          Real tlsic = hugn/2*Usupg.size()*zeta;

          //Calculate CFL
          Real cfl = fabs(u*dt/elem->hmin())+fabs(v*dt/elem->hmin());
          if (maxcfl<cfl)
              maxcfl = cfl;
          if (maxtsupg<tsupg)
              maxtsupg = tsupg;

         //tsupg = 0;
         // forcex = 0;
         // forcey = 0;


          if(flag==1)//Residual calculation
          {

            for (unsigned int i=0;i<n_u_dofs;i++)
            {
            Reu(i) += JxW[qp]*(
                                +(Rho*(u - uold)*phi[i][qp])/dt
                                +(2*dudx*(dphi[i][qp](0))*Mu)/Re + ((dudy + dvdx)*(dphi[i][qp](1))*Mu)/Re
                                +Rho*(dudx*u + dudy*v)*phi[i][qp]
                                +-(dphi[i][qp](0))*p
                                +0
                                +-(((1./We)*kurvature*grad_H(0))*phi[i][qp])
                                //+-(1./We)*(-jump*dphi[i][qp](0) + forcex*phi[i][qp] )
                                +(Rho*tsupg*(u - uold)*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/dt
                                +-(((d2udx2 + d2udy2)*Mu*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/Re)
                                +Rho*tsupg*(dudx*u + dudy*v)*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                                +dpdx*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                                +-(((1./We)*kurvature*grad_H(0))*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))
                                //+-((1./We)*(grad_jump(0)+forcex)*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))
                               );

            Rev(i) += JxW[qp]*(
                                +(Rho*(v - vold)*phi[i][qp])/dt
                                +((dudy + dvdx)*(dphi[i][qp](0))*Mu)/Re + (2*dvdy*(dphi[i][qp](1))*Mu)/Re
                                +Rho*(dvdx*u + dvdy*v)*phi[i][qp]
                                +-(dphi[i][qp](1))*p
                                +0
                                +-(((1./We)*kurvature*grad_H(1) + Rho*(-1./(Fr*Fr))*gravity_on)*phi[i][qp])
                                //+-(1./We)*(-jump*dphi[i][qp](1) + forcey*phi[i][qp] ) + -((Rho*(-1./(Fr*Fr))*gravity_on)*phi[i][qp])
                                +(Rho*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)*(v - vold))/dt
                                +-(((d2vdx2 + d2vdy2)*Mu*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))/Re)
                                +Rho*tsupg*(dvdx*u + dvdy*v)*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                                +dpdy*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v)
                                +-(((1./We)*kurvature*grad_H(1) + Rho*(-1./(Fr*Fr))*gravity_on)*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))
                                //+-(((1./We)*(grad_jump(1)+forcey) + Rho*(-1./(Fr*Fr))*gravity_on)*tsupg*((dphi[i][qp](0))*u + (dphi[i][qp](1))*v))
                               );

            }

            for (unsigned int i=0;i<n_p_dofs;i++)
            {
            Rep(i) += JxW[qp]*(
                                +0
                                +0
                                +0
                                +0
                                +(dudx + dvdy)*psi[i][qp]
                                +0
                                +0
                                +0
                                +0
                                +0
                                +0
                               );

            }



          }//end if flag==1 (residual calculation



          if(flag==2)
          {



              for (unsigned int i = 0; i < n_u_dofs; i++) {

                for (unsigned int j = 0; j < n_u_dofs; j++) {
                    Kuu(i, j) += JxW[qp]*(
                            +(Rho * phi[j][qp] * phi[i][qp]) / dt
                            + ((2 * (dphi[j][qp](0))*(dphi[i][qp](0)) + (dphi[j][qp](1))*(dphi[i][qp](1))) * Mu) / Re
                            + Rho * ((dphi[j][qp](0)) * u + dudx * phi[j][qp] + (dphi[j][qp](1)) * v) * phi[i][qp]
                            + 0
                            + 0
                            + 0
                            + (Rho * tsupg * phi[j][qp]*(2 * (dphi[i][qp](0)) * u - (dphi[i][qp](0)) * uold + (dphi[i][qp](1)) * v)) / dt
                            + -((Mu * tsupg * ((d2phidx2[j][qp] + d2phidy2[j][qp])*(dphi[i][qp](0)) * u + (d2udx2 + d2udy2)*(dphi[i][qp](0)) * phi[j][qp] +
                            (d2phidx2[j][qp] + d2phidy2[j][qp])*(dphi[i][qp](1)) * v)) / Re)
                            + Rho * tsupg * ((dphi[j][qp](0)) * u * ((dphi[i][qp](0)) * u + (dphi[i][qp](1)) * v) + dudx * phi[j][qp]*(2 * (dphi[i][qp](0)) * u + (dphi[i][qp](1)) * v) +
                            v * ((dphi[j][qp](1))*(dphi[i][qp](0)) * u + dudy * (dphi[i][qp](0)) * phi[j][qp] + (dphi[j][qp](1))*(dphi[i][qp](1)) * v))
                            + dpdx * (dphi[i][qp](0)) * tsupg * phi[j][qp]
                            +-((dphi[i][qp](0))*((1. / We) * kurvature * grad_H(0)) * tsupg * phi[j][qp])
                            //+-((dphi[i][qp](0))*((1. / We) * (grad_jump(0)+forcex)) * tsupg * phi[j][qp])
                            );

                    Kuv(i, j) += JxW[qp]*(
                            +0
                            + ((dphi[j][qp](0))*(dphi[i][qp](1)) * Mu) / Re
                            + dudy * Rho * phi[j][qp] * phi[i][qp]
                            + 0
                            + 0
                            + 0
                            + ((dphi[i][qp](1)) * Rho * tsupg * (u - uold) * phi[j][qp]) / dt
                            + -(((d2udx2 + d2udy2)*(dphi[i][qp](1)) * Mu * tsupg * phi[j][qp]) / Re)
                            + Rho * tsupg * (dudy * (dphi[i][qp](0)) * u + dudx * (dphi[i][qp](1)) * u + 2 * dudy * (dphi[i][qp](1)) * v) * phi[j][qp]
                            + dpdx * (dphi[i][qp](1)) * tsupg * phi[j][qp]
                            +-((dphi[i][qp](1))*((1. / We) * (kurvature) * grad_H(0)) * tsupg * phi[j][qp])
                            //+-((dphi[i][qp](1))*((1. / We) * (grad_jump(0)+forcex)) * tsupg * phi[j][qp])
                            );

                    Kvu(i, j) += JxW[qp]*(
                            +0
                            + ((dphi[j][qp](1))*(dphi[i][qp](0)) * Mu) / Re
                            + dvdx * Rho * phi[j][qp] * phi[i][qp]
                            + 0
                            + 0
                            + 0
                            + ((dphi[i][qp](0)) * Rho * tsupg * phi[j][qp]*(v - vold)) / dt
                            + -(((d2vdx2 + d2vdy2)*(dphi[i][qp](0)) * Mu * tsupg * phi[j][qp]) / Re)
                            + Rho * tsupg * phi[j][qp]*(2 * dvdx * (dphi[i][qp](0)) * u + dvdy * (dphi[i][qp](0)) * v + dvdx * (dphi[i][qp](1)) * v)
                            + dpdy * (dphi[i][qp](0)) * tsupg * phi[j][qp]
                            //+-((dphi[i][qp](0))*((1. / We) *(kurvature) * grad_H(1) + Rho * (-1. / (Fr * Fr)) * gravity_on) * tsupg * phi[j][qp])
                            +-((dphi[i][qp](0))*((1. / We) *(grad_jump(1)+forcey) + Rho * (-1. / (Fr * Fr)) * gravity_on) * tsupg * phi[j][qp])
                            );

                    Kvv(i, j) += JxW[qp]*(
                            +(Rho * phi[j][qp] * phi[i][qp]) / dt
                            + (((dphi[j][qp](0))*(dphi[i][qp](0)) + 2 * (dphi[j][qp](1))*(dphi[i][qp](1))) * Mu) / Re
                            + Rho * ((dphi[j][qp](0)) * u + (dphi[j][qp](1)) * v + dvdy * phi[j][qp]) * phi[i][qp]
                            + 0
                            + 0
                            + 0
                            + (Rho * tsupg * phi[j][qp]*((dphi[i][qp](0)) * u + 2 * (dphi[i][qp](1)) * v - (dphi[i][qp](1)) * vold)) / dt
                            + -((Mu * tsupg * ((d2phidx2[j][qp] + d2phidy2[j][qp])*((dphi[i][qp](0)) * u + (dphi[i][qp](1)) * v) +
                            (d2vdx2 + d2vdy2)*(dphi[i][qp](1)) * phi[j][qp])) / Re)
                            + Rho * tsupg * (((dphi[j][qp](0)) * u + (dphi[j][qp](1)) * v)*((dphi[i][qp](0)) * u + (dphi[i][qp](1)) * v) +
                            (dvdy * (dphi[i][qp](0)) * u + dvdx * (dphi[i][qp](1)) * u + 2 * dvdy * (dphi[i][qp](1)) * v) * phi[j][qp])
                            + dpdy * (dphi[i][qp](1)) * tsupg * phi[j][qp]
                            +-((dphi[i][qp](1))*((1. / We) *(kurvature) * grad_H(1) + Rho * (-1. / (Fr * Fr)) * gravity_on) * tsupg * phi[j][qp])
                            //+-((dphi[i][qp](1))*((1. / We) *(grad_jump(1)+forcey) + Rho * (-1. / (Fr * Fr)) * gravity_on) * tsupg * phi[j][qp])
                            );

                }

                for (unsigned int j = 0; j < n_p_dofs; j++) {
                    Kup(i, j) += JxW[qp]*(
                            +0
                            + 0
                            + 0
                            + -(dphi[i][qp](0)) * psi[j][qp]
                            + 0
                            + 0
                            + 0
                            + 0
                            + 0
                            + (dpsi[j][qp](0)) * tsupg * ((dphi[i][qp](0)) * u + (dphi[i][qp](1)) * v)
                            + 0
                            );

                    Kvp(i, j) += JxW[qp]*(
                            +0
                            + 0
                            + 0
                            + -(dphi[i][qp](1)) * psi[j][qp]
                            + 0
                            + 0
                            + 0
                            + 0
                            + 0
                            + (dpsi[j][qp](1)) * tsupg * ((dphi[i][qp](0)) * u + (dphi[i][qp](1)) * v)
                            + 0
                            );

                }
            }

            for (unsigned int i = 0; i < n_p_dofs; i++) {
                for (unsigned int j = 0; j < n_u_dofs; j++) {
                    Kpu(i, j) += JxW[qp]*(
                            +0
                            + 0
                            + 0
                            + 0
                            + (dphi[j][qp](0)) * psi[i][qp]
                            + 0
                            + 0
                            + 0
                            + 0
                            + 0
                            + 0
                            );

                    Kpv(i, j) += JxW[qp]*(
                            +0
                            + 0
                            + 0
                            + 0
                            + (dphi[j][qp](1)) * psi[i][qp]
                            + 0
                            + 0
                            + 0
                            + 0
                            + 0
                            + 0
                            );

                }
            }





          }//end if flag==2 (jacobian calculation)


      }//End of quadrature loop

       //TODO Boundaries genericos
      RealVectorValue (*stokes_dv_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< RealVectorValue (*)(const Real , const Real , const Real)>("stokes_dv_boundary_fptr");
      Real (*stokes_dp_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< Real (*)(const Real , const Real , const Real)>("stokes_dp_boundary_fptr");
      bool (*stokes_is_dv_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< bool (*)(const Real , const Real , const Real)>("stokes_is_dv_boundary_fptr");
      bool (*stokes_is_dp_boundary_fptr)(const Real , const Real , const Real) =
              es.parameters.get< bool (*)(const Real , const Real , const Real)>("stokes_is_dp_boundary_fptr");

      //Begin boundary conditions
      if (flag==1)//Residual boundary conditions
      {

        const Real penalty = 1.e10;
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == NULL)
            {

              AutoPtr<Elem> side (elem->build_side(s));

              // Loop over the nodes on the side.
              for (unsigned int ns=0; ns< side->n_nodes(); ns++)
                {

                  // Loop over the nodes of the element
                  for (unsigned int n=0; n<elem->n_nodes(); n++)
                  {
                    if (elem->node(n) == side->node(ns))
                    {

                        Node* sidenode = elem->get_node(n);
                        Real x = sidenode->operator ()(0);
                        Real y = sidenode->operator ()(1);
                        Real z = sidenode->operator ()(2);

                        Real u = soln(dof_indices_u[n]);
                        Real v = soln(dof_indices_v[n]);
                        //TODO 3D
                        Real p = soln(dof_indices_p[n]);


                       if (stokes_dv_boundary_fptr!=NULL && stokes_is_dv_boundary_fptr!=NULL)
                       {

                            if (n<n_u_dofs && (*stokes_is_dv_boundary_fptr)(x,y,z)==true)
                            {
                                const RealVectorValue vector = (*stokes_dv_boundary_fptr)(x,y,z);

                                Reu(n) = penalty*(u - vector(0));
                                Rev(n) = penalty*(v - vector(1));

                                //TODO 3D
                            }
                       }
                       else
                       {
                               if (n<n_u_dofs)
                               {

                                Reu(n) = penalty*(u);
                                Rev(n) = penalty*(v);

                                //TODO 3D
                               }

                       }

                       if (stokes_dp_boundary_fptr!=NULL && stokes_is_dp_boundary_fptr!=NULL)
                       {
                            if (n<n_p_dofs && (*stokes_is_dp_boundary_fptr)(x,y,z)==true)
                            {
                                const Real pressure = (*stokes_dp_boundary_fptr)(x,y,z);

                                Rep(n) = penalty*(p - pressure);

                            }
                       }
                       else
                       {
                            const unsigned int pressure_node = 0;
                            if (elem->node(n) == pressure_node)
                            {
                              const Real p_value = 0.0;
                                Rep(n) = penalty*(p - p_value);
                            }

                       }

                    }

                  }
                } // end face node loop
            } // end if (elem->neighbor(side) == NULL)

      }// end residual boundary condition section



      //Begin jacobian boundary conditions
      if(flag==2)
      {

        const Real penalty = 1.e10;
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == NULL)
            {

              AutoPtr<Elem> side (elem->build_side(s));

              // Loop over the nodes on the side.
              for (unsigned int ns=0; ns< side->n_nodes(); ns++)
                {

                  // Loop over the nodes of the element
                  for (unsigned int n=0; n<elem->n_nodes(); n++)
                  {
                    if (elem->node(n) == side->node(ns))
                    {

                        Node* sidenode = elem->get_node(n);
                        Real x = sidenode->operator ()(0);
                        Real y = sidenode->operator ()(1);
                        Real z = sidenode->operator ()(2);

                       if (stokes_dv_boundary_fptr!=NULL && stokes_is_dv_boundary_fptr!=NULL)
                       {

                            if (n<n_u_dofs && (*stokes_is_dv_boundary_fptr)(x,y,z)==true)
                            {
                                const RealVectorValue vector = (*stokes_dv_boundary_fptr)(x,y,z);
                                // Matrix contribution.
                                Kuu(n,n) += penalty;
                                Kvv(n,n) += penalty;

                                //TODO 3D
                            }
                       }
                       else
                       {
                               if (n<n_u_dofs)
                               {
                                Kuu(n,n) += penalty;
                                Kvv(n,n) += penalty;

                                //TODO 3D
                               }

                       }

                       if (stokes_dp_boundary_fptr!=NULL && stokes_is_dp_boundary_fptr!=NULL)
                       {
                            if (n<n_p_dofs && (*stokes_is_dp_boundary_fptr)(x,y,z)==true)
                            {
                                const Real pressure = (*stokes_dp_boundary_fptr)(x,y,z);

                                Kpp(n,n) += penalty;

                            }
                       }
                       else
                       {
                            const unsigned int pressure_node = 0;
                            if (elem->node(n) == pressure_node)
                            {
                              const Real p_value = 0.0;
                              Kpp(n,n) += penalty;
                            }

                       }

                    }

                  }
                } // end face node loop
            } // end if (elem->neighbor(side) == NULL)

      }
      // end jacobian boundary condition section (flag==2)


      if(flag==1)
      {
      dof_map.constrain_element_vector (Residual, dof_indices);
      residual.add_vector (Residual, dof_indices);
      }

      if(flag==2)
      {
      dof_map.constrain_element_matrix (Ke, dof_indices);
      jacobian.add_matrix (Ke, dof_indices);
      }

    }//end of element loop

  es.parameters.set<Real>("CFL")=maxcfl;
  es.parameters.set<Real>("maxtsupg")=maxtsupg;

}



void assemble_jump (EquationSystems& es, //TODO functions is deprecated
                      const std::string& system_name)
{

  libmesh_assert (system_name == "Jump Pres System"
                  || system_name == "Jump Vel System");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_jump =
          es.get_system<LinearImplicitSystem> (system_name);

    LinearImplicitSystem & system_wforce =
          es.get_system<LinearImplicitSystem> ("Whole Force System");

  const unsigned int wforcex_var = system_wforce.variable_number ("wforcex");
  const unsigned int wforcey_var = system_wforce.variable_number ("wforcey");
  //const unsigned int forcez_var = system_force.variable_number ("forcez");
  const unsigned int wforcez_var = (dim==3) ? system_wforce.variable_number ("wforcez") : 0;

  const DofMap & dof_map = system_jump.get_dof_map();
  const DofMap & dof_map_wforce = system_wforce.get_dof_map();

  FEType fe_pres_type = system_jump.variable_type(0);
  FEType fe_vel_type = system_wforce.variable_type(0);

  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  QGauss qrule (dim, fe_pres_type.default_quadrature_order());
  fe_pres->attach_quadrature_rule (&qrule);
  fe_vel->attach_quadrature_rule (&qrule);
  
  AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));
  QGauss qface (dim-1,fe_pres_type.default_quadrature_order());
  fe_pres_face->attach_quadrature_rule (&qface);



  const std::vector<Real>& JxW = fe_pres->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();


    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;


  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_wforcex;
  std::vector<unsigned int> dof_indices_wforcey;
  std::vector<unsigned int> dof_indices_wforcez;

  //Interface * iop = es.parameters.get<Interface*>("interface_obj_ptr");

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el ; ++el)
      {
        const Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices);

        fe_pres->reinit (elem);
        fe_vel->reinit(elem);

          dof_map_wforce.dof_indices (elem, dof_indices_wforcex, wforcex_var);
          dof_map_wforce.dof_indices (elem, dof_indices_wforcey, wforcey_var);
          if (dim==3)
              dof_map_wforce.dof_indices (elem, dof_indices_wforcez, wforcez_var);


        Ke.resize (dof_indices.size(),
                   dof_indices.size());

        Fe.resize (dof_indices.size());


      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {

          Number wforcex = 0;
          Number wforcey = 0;
          Number wforcez = 0;

          for (unsigned int l=0; l<phi.size(); l++)
          {
              wforcex += phi[l][qp]*system_wforce.current_solution(dof_indices_wforcex[l]);
              wforcey += phi[l][qp]*system_wforce.current_solution(dof_indices_wforcey[l]);
              if(dim==3)
                  wforcez += phi[l][qp]*system_wforce.current_solution(dof_indices_wforcez[l]);

          }
          
          RealVectorValue wforce(wforcex,wforcey,wforcez);


          for (unsigned int i=0; i<psi.size(); i++)
          {
              // The RHS contribution
              //Fe(i) += -JxW[qp]*(-kurvature*dH*(-normal*dpsi[i][qp]));
              Fe(i) += JxW[qp]*((wforce*dpsi[i][qp]));

              for (unsigned int j=0; j<psi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(dpsi[i][qp]*dpsi[j][qp] );
                }

          }
      }//end or quadrature points
      

        {//Boundary section

          for (unsigned int side=0; side<elem->n_sides(); side++)
            if (elem->neighbor(side) == NULL)
              {
                const std::vector<std::vector<Real> >&  psi_face = fe_pres_face->get_phi();
                const std::vector<Real>& JxW_face = fe_pres_face->get_JxW();

                fe_pres_face->reinit(elem, side);

                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                  {

                    const Real penalty = 1.e10;
                    const Real value = 0;

                    for (unsigned int i=0; i<psi_face.size(); i++)
                      for (unsigned int j=0; j<psi_face.size(); j++)
                        Ke(i,j) += JxW_face[qp]*penalty*psi_face[i][qp]*psi_face[j][qp];

                    for (unsigned int i=0; i<psi_face.size(); i++)
                      Fe(i) += JxW_face[qp]*penalty*value*psi_face[i][qp];
                  }
              }
        }


        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        system_jump.matrix->add_matrix (Ke, dof_indices);
        system_jump.rhs->add_vector    (Fe, dof_indices);


    } // end of element loop

  // That's it.
  return;
}

void assemble_force (EquationSystems& es, //TODO functions is deprecated
                      const std::string& system_name)
{

  libmesh_assert (system_name == "Force System");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_force =
          es.get_system<LinearImplicitSystem> ("Force System");
  LinearImplicitSystem & system_wforce =
          es.get_system<LinearImplicitSystem> ("Whole Force System");
  LinearImplicitSystem & system_jump =
          es.get_system<LinearImplicitSystem> ("Jump Vel System");
  TransientAdvectionDiffusionSystem & system_normalx =
    es.get_system<TransientAdvectionDiffusionSystem> ("normalx SYSTEM");
  TransientAdvectionDiffusionSystem & system_normaly =
    es.get_system<TransientAdvectionDiffusionSystem> ("normaly SYSTEM");
  TransientAdvectionDiffusionSystem & system_normalz =
    es.get_system<TransientAdvectionDiffusionSystem> ("normalz SYSTEM");

  LinearImplicitSystem & system_heaviside =
    es.get_system<LinearImplicitSystem> ("Heaviside Function");

  
  const unsigned int forcex_var = system_force.variable_number ("forcex");
  const unsigned int forcey_var = system_force.variable_number ("forcey");
  const unsigned int forcez_var = (dim==3) ? system_force.variable_number ("forcez") : 0;

  const unsigned int wforcex_var = system_wforce.variable_number ("wforcex");
  const unsigned int wforcey_var = system_wforce.variable_number ("wforcey");
  const unsigned int wforcez_var = (dim==3) ? system_wforce.variable_number ("wforcez") : 0;



  FEType fe_vel_type = system_force.variable_type(0);
  FEType fe_pres_type = system_heaviside.variable_type(0);

  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  fe_pres->attach_quadrature_rule (&qrule);
  fe_vel->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe_vel->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();


  const DofMap & dof_map = system_force.get_dof_map();
  const DofMap & dof_map_wforce = system_wforce.get_dof_map();
  const DofMap & dof_map_jump = system_jump.get_dof_map();
  const DofMap & dof_map_heaviside = system_heaviside.get_dof_map();
  const DofMap & dof_map_normalx = system_normalx.get_dof_map();
  const DofMap & dof_map_normaly = system_normaly.get_dof_map();
  const DofMap & dof_map_normalz = system_normalz.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke),Kuw(Ke),
    Kvu(Ke), Kvv(Ke),Kvw(Ke),
    Kwu(Ke), Kwv(Ke),Kww(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_forcex;
  std::vector<unsigned int> dof_indices_forcey;
  std::vector<unsigned int> dof_indices_forcez;
  std::vector<unsigned int> dof_indices_wforcex;
  std::vector<unsigned int> dof_indices_wforcey;
  std::vector<unsigned int> dof_indices_wforcez;
  std::vector<unsigned int> dof_indices_normalx;
  std::vector<unsigned int> dof_indices_normaly;
  std::vector<unsigned int> dof_indices_normalz;
  std::vector<unsigned int> dof_indices_jump;
  std::vector<unsigned int> dof_indices_heaviside;

  Interface * iop = es.parameters.get<Interface*>("interface_obj_ptr");

  Real time = es.parameters.get<Real>("time");
  Real dt = es.parameters.get<Real>("dt");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  Number RhoA =es.parameters.get<Real>("RhoA"); //Density
  Number RhoB =es.parameters.get<Real>("RhoB");; //Density


  for ( ; el != end_el; ++el)
  {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_forcex, forcex_var);
      dof_map.dof_indices (elem, dof_indices_forcey, forcey_var);
      if (dim==3)
          dof_map.dof_indices (elem, dof_indices_forcez, forcez_var);

      dof_map_wforce.dof_indices (elem, dof_indices_wforcex, wforcex_var);
      dof_map_wforce.dof_indices (elem, dof_indices_wforcey, wforcey_var);
      if (dim==3)
          dof_map_wforce.dof_indices (elem, dof_indices_wforcez, wforcez_var);

      dof_map_jump.dof_indices(elem,dof_indices_jump);
      dof_map_heaviside.dof_indices(elem,dof_indices_heaviside);
      dof_map_normalx.dof_indices(elem,dof_indices_normalx);
      dof_map_normaly.dof_indices(elem,dof_indices_normaly);
      if (dim==3)
          dof_map_normalz.dof_indices(elem,dof_indices_normalz);


      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_forcex.size();
      const unsigned int n_heaviside_dofs = dof_indices_heaviside.size();

      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (forcex_var*n_u_dofs, forcex_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (forcex_var*n_u_dofs, forcey_var*n_u_dofs, n_u_dofs, n_u_dofs);
      if (dim==3)
          Kuw.reposition (forcex_var*n_u_dofs, forcez_var*n_u_dofs, n_u_dofs, n_u_dofs);

      Kvu.reposition (forcey_var*n_u_dofs, forcex_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kvv.reposition (forcey_var*n_u_dofs, forcey_var*n_u_dofs, n_u_dofs, n_u_dofs);
      if (dim==3)
          Kvw.reposition (forcey_var*n_u_dofs, forcez_var*n_u_dofs, n_u_dofs, n_u_dofs);

      if (dim==3)
      {
          Kwu.reposition (forcez_var*n_u_dofs, forcex_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kwv.reposition (forcez_var*n_u_dofs, forcey_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kww.reposition (forcez_var*n_u_dofs, forcez_var*n_u_dofs, n_u_dofs, n_u_dofs);
      }


      Fu.reposition (forcex_var*n_u_dofs, n_u_dofs);
      Fv.reposition (forcey_var*n_u_dofs, n_u_dofs);
      if (dim==3)
          Fw.reposition (forcez_var*n_u_dofs, n_u_dofs);


      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          Number wforcex = 0;
          Number wforcey = 0;
          Number wforcez = 0;
          Number normalx = 0;
          Number normaly = 0;
          Number normalz = 0;
          Number jump = 0;
          Gradient grad_jump;
          Gradient grad_H;
          Number heaviside =0;



          for (unsigned int l=0; l<phi.size(); l++)
          {

              wforcex += phi[l][qp]*system_wforce.current_solution(dof_indices_wforcex[l]);
              wforcey += phi[l][qp]*system_wforce.current_solution(dof_indices_wforcey[l]);
              if(dim==3)
                  wforcez += phi[l][qp]*system_wforce.current_solution(dof_indices_wforcez[l]);

              normalx += phi[l][qp]*system_normalx.current_solution(dof_indices_normalx[l]);
              normaly += phi[l][qp]*system_normaly.current_solution(dof_indices_normaly[l]);
              if(dim==3)
                  normalz += phi[l][qp]*system_normalz.current_solution(dof_indices_normalz[l]);

          }

          for (unsigned int l=0; l<phi.size(); l++)
          {
              grad_jump.add_scaled (dphi[l][qp],system_jump.current_solution (dof_indices_jump[l]));
              //jump += psi[l][qp]*system_jump.current_solution(dof_indices_jump[l]);
          }

          grad_H.zero();
          for(unsigned int l=0; l<n_heaviside_dofs; l++)
            {
                heaviside += psi[l][qp]*system_heaviside.current_solution(dof_indices_heaviside[l]);
                grad_H.add_scaled (dpsi[l][qp],system_heaviside.current_solution (dof_indices_heaviside[l]));
            }


          Real diff=0.0001;
          Real tsupg = 0.0005;
          Real Rho = RhoB + (RhoA - RhoB)*heaviside;
          RealVectorValue seminormal (normalx,normaly,normalz);
          seminormal.unit();




          for (unsigned int i=0; i<phi.size(); i++)
          {
              // The RHS contribution

              //Fu(i) += JxW[qp]*(-kurvature*dH*normal(0)-grad_jump(0))*phi[i][qp]//;//-grad_jump(0)
              Fu(i) += JxW[qp]*(wforcex - grad_jump(0))*phi[i][qp];//-grad_jump(0)

              //Fv(i) += JxW[qp]*(-kurvature*dH*normal(1) -grad_jump(1))*phi[i][qp];//-grad_jump(1)
              Fv(i) += JxW[qp]*(wforcey - grad_jump(1))*phi[i][qp];//-grad_jump(1)

              if (dim==3)
                  Fw(i) += JxW[qp]*(wforcez - grad_jump(2))*phi[i][qp];//-grad_jump(2)

              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            +diff*dphi[i][qp]*dphi[j][qp]
                                            +tsupg*(seminormal*dphi[i][qp])*(seminormal*dphi[j][qp])
                                        );
                  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            +diff*dphi[i][qp]*dphi[j][qp]
                                            +tsupg*(seminormal*dphi[i][qp])*(seminormal*dphi[j][qp])
                                        );
                  if(dim==3)
                      Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            +diff*dphi[i][qp]*dphi[j][qp]
                                            +tsupg*(seminormal*dphi[i][qp])*(seminormal*dphi[j][qp])
                                        );

                }

          }


      }//end or quadrature points


        //BEGIN BOUNDARY CONDITIONS SECTION
        /*{
          const Real penalty = 1.e10;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            if (elem->neighbor(s) == NULL)
              {

                AutoPtr<Elem> side (elem->build_side(s));

                for (unsigned int ns=0; ns<side->n_nodes(); ns++)
                  {


                    for (unsigned int n=0; n<elem->n_nodes(); n++)
                      if (elem->node(n) == side->node(ns))
                        {
                          Kuu(n,n) += penalty;
                          Kvv(n,n) += penalty;
                          Kww(n,n) += penalty;

                          Fu(n) += penalty*0;
                          Fv(n) += penalty*0;
                          Fw(n) += penalty*0;
                        }
                  } // end face node loop
              } // end if (elem->neighbor(side) == NULL)

        } // end boundary condition section*/



      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system_force.matrix->add_matrix (Ke, dof_indices);
      system_force.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

  // That's it.
  return;
}

void assemble_whole_force (EquationSystems& es, //TODO functions is deprecated
                      const std::string& system_name)
{

  libmesh_assert (system_name == "Whole Force System");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system_wforce =
          es.get_system<LinearImplicitSystem> ("Whole Force System");

  LinearImplicitSystem & system_kurvature =
    es.get_system<LinearImplicitSystem> ("Kurvature");
  TransientAdvectionDiffusionSystem & system_normalx =
    es.get_system<TransientAdvectionDiffusionSystem> ("normalx SYSTEM");
  TransientAdvectionDiffusionSystem & system_normaly =
    es.get_system<TransientAdvectionDiffusionSystem> ("normaly SYSTEM");
  TransientAdvectionDiffusionSystem & system_normalz =
    es.get_system<TransientAdvectionDiffusionSystem> ("normalz SYSTEM");
  TransientLinearImplicitSystem & system_reinit =
    es.get_system<TransientLinearImplicitSystem> ("Reinit_LevelSet");
  LinearImplicitSystem & system_heaviside =
    es.get_system<LinearImplicitSystem> ("Heaviside Function");



  // Numeric ids corresponding to each variable in the system
  const unsigned int forcex_var = system_wforce.variable_number ("wforcex");
  const unsigned int forcey_var = system_wforce.variable_number ("wforcey");
  const unsigned int forcez_var = (dim==3) ? system_wforce.variable_number ("wforcez") : 0;


  FEType fe_vel_type = system_wforce.variable_type(0);

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());


  fe_vel->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe_vel->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();


  const DofMap & dof_map = system_wforce.get_dof_map();
  const DofMap & dof_map_normalx = system_normalx.get_dof_map();
  const DofMap & dof_map_normaly = system_normaly.get_dof_map();
  const DofMap & dof_map_normalz = system_normalz.get_dof_map();
  const DofMap & dof_map_kurvature = system_kurvature.get_dof_map();
  const DofMap & dof_map_reinit = system_reinit.get_dof_map();
  const DofMap & dof_map_heaviside = system_heaviside.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke),Kuw(Ke),
    Kvu(Ke), Kvv(Ke),Kvw(Ke),
    Kwu(Ke), Kwv(Ke),Kww(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_forcex;
  std::vector<unsigned int> dof_indices_forcey;
  std::vector<unsigned int> dof_indices_forcez;
  std::vector<unsigned int> dof_indices_normalx;
  std::vector<unsigned int> dof_indices_normaly;
  std::vector<unsigned int> dof_indices_normalz;
  std::vector<unsigned int> dof_indices_kurvature;
  std::vector<unsigned int> dof_indices_reinit;
  std::vector<unsigned int> dof_indices_heaviside;

  Interface * iop = es.parameters.get<Interface*>("interface_obj_ptr");

  Real time = es.parameters.get<Real>("time");
  Real dt = es.parameters.get<Real>("dt");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_forcex, forcex_var);
      dof_map.dof_indices (elem, dof_indices_forcey, forcey_var);
      if (dim==3)
          dof_map.dof_indices (elem, dof_indices_forcez, forcez_var);

      dof_map_normalx.dof_indices (elem, dof_indices_normalx);
      dof_map_normaly.dof_indices (elem, dof_indices_normaly);
      if (dim==3)
          dof_map_normalz.dof_indices (elem, dof_indices_normalz);
      dof_map_reinit.dof_indices (elem, dof_indices_reinit);
      dof_map_kurvature.dof_indices (elem, dof_indices_kurvature);
      dof_map_heaviside.dof_indices (elem, dof_indices_heaviside);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_forcex.size();

      fe_vel->reinit  (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (forcex_var*n_u_dofs, forcex_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (forcex_var*n_u_dofs, forcey_var*n_u_dofs, n_u_dofs, n_u_dofs);
      if (dim==3)
          Kuw.reposition (forcex_var*n_u_dofs, forcez_var*n_u_dofs, n_u_dofs, n_u_dofs);

      Kvu.reposition (forcey_var*n_u_dofs, forcex_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kvv.reposition (forcey_var*n_u_dofs, forcey_var*n_u_dofs, n_u_dofs, n_u_dofs);
      if (dim==3)
          Kvw.reposition (forcey_var*n_u_dofs, forcez_var*n_u_dofs, n_u_dofs, n_u_dofs);

      if (dim==3)
      {
          Kwu.reposition (forcez_var*n_u_dofs, forcex_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kwv.reposition (forcez_var*n_u_dofs, forcey_var*n_u_dofs, n_u_dofs, n_u_dofs);
          Kww.reposition (forcez_var*n_u_dofs, forcez_var*n_u_dofs, n_u_dofs, n_u_dofs);
      }

      Fu.reposition (forcex_var*n_u_dofs, n_u_dofs);
      Fv.reposition (forcey_var*n_u_dofs, n_u_dofs);
      if (dim==3)
          Fw.reposition (forcez_var*n_u_dofs, n_u_dofs);


      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Values to hold the solution & its gradient at the previous timestep.
          Number lvlset_value = 0;
          Number kurvature = 0;
          Number normalx = 0;
          Number normaly = 0;
          Number normalz = 0;
          Number heaviside = 0;
          Gradient grad_H;

          for (unsigned int l=0; l<phi.size(); l++)
          {
              lvlset_value += phi[l][qp]*system_reinit.current_solution(dof_indices_reinit[l]);
              kurvature += phi[l][qp]*system_kurvature.current_solution(dof_indices_kurvature[l]);
              normalx += phi[l][qp]*system_normalx.current_solution(dof_indices_normalx[l]);
              normaly += phi[l][qp]*system_normaly.current_solution(dof_indices_normaly[l]);
              if(dim==3)
                  normalz += phi[l][qp]*system_normalz.current_solution(dof_indices_normalz[l]);

          }

          for(unsigned int l=0; l<dof_indices_heaviside.size(); l++)
            {
                heaviside += phi[l][qp]*system_heaviside.current_solution(dof_indices_heaviside[l]);
                grad_H.add_scaled (dphi[l][qp],system_heaviside.current_solution (dof_indices_heaviside[l]));
            }

          RealVectorValue normal;
          normal(0)=normalx;
          normal(1)=normaly;
          if (dim==3)
              normal(2)=normalz;
          normal.unit();
          Real dH = iop->dH(lvlset_value);
          Real epsilon = iop->getEpsilon();
          Real diff= 0.0;//0.0004;
          Real diffdH=0.0;


          for (unsigned int i=0; i<phi.size(); i++)
          {
              // The RHS contribution

              //Fu(i) += JxW[qp]*(kurvature*dH*normal(0)*phi[i][qp]);//-grad_jump(0)
              Fu(i) += JxW[qp]*(kurvature*grad_H(0)*phi[i][qp]);//-grad_jump(0)
              

              //Fv(i) += JxW[qp]*(kurvature*dH*normal(1))*phi[i][qp];//-grad_jump(1)
              Fv(i) += JxW[qp]*(kurvature*grad_H(1)*phi[i][qp]);//-grad_jump(1)

              if (dim==3)
                  Fw(i) += JxW[qp]*(kurvature*dH*normal(2))*phi[i][qp];//-grad_jump(2)

              for (unsigned int j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            +diff*dphi[i][qp]*dphi[j][qp]*(1-dH*epsilon)
                                            +diffdH*dphi[i][qp]*dphi[j][qp]*dH
                                        );
                  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            +diff*dphi[i][qp]*dphi[j][qp]*(1-dH*epsilon)
                                            +diffdH*dphi[i][qp]*dphi[j][qp]*dH
                                        );
                  if(dim==3)
                      Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]
                                            +diff*dphi[i][qp]*dphi[j][qp]
                                            +diffdH*dphi[i][qp]*dphi[j][qp]*dH
                                        );

                }

          }


      }//end or quadrature points



      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system_wforce.matrix->add_matrix (Ke, dof_indices);
      system_wforce.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

  // That's it.
  return;
}

