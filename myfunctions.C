
#include "myfunctions.h"


Number function_eliptic_bubble (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
    Real x = p(0);
    Real y = p(1);
    Real num2=0;
    Real Rx = 0.5;
    Real Ry = 0.525;
    Real x0 = 0;
    Real y0 = 0;
    Real alpha = sqrt(pow((x-x0)/Rx,2)+pow((y-y0)/Ry,2));
    Real theta;
    if (x==0 && y==0)
        theta = -asin(1);
    else
        theta = atan(x/y);

    num2 = (1-alpha)*sqrt(pow(Rx*cos(theta),2)+pow(Ry*sin(theta),2));
    return num2;
}

Number function_two_bubbles (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
    Real x = p(0);
    Real y = p(1);
    Real z = 0;
    Real num2=0;
    Real R = 0.4;
    Real x0 = 0;
    Real y0 = 1;
    Real r = sqrt(pow((x-x0),2)+pow((y-y0),2));
    if (y<=1.45)
        num2 = (R - r)*1;
    //second bubble
    R=0.5;
    x0=0;
    y0=2;
    r = sqrt(pow((x-x0),2)+pow((y-y0),2));
    if (y>1.45)
        num2 = (R - r)*1;
    //surface
    Real y_limit = 3;
    Real alpha = (y-2)/(y_limit-2);
    if (y>2.7)
        num2 = (y - y_limit)*alpha + (R-r)*(1-alpha);

    if (y>y_limit)
        num2 = y - y_limit;


  //return exp(-num/den)/(4.*t + 1.);
  return num2;

}
Number function_one_bubble (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
    Real x = p(0);
    Real y = p(1);
    Real z = 0;
    Real num2=0;
    Real R = 0.5;
    Real x0 = 0.0;
    Real y0 = 0.0;
    Real r = sqrt(pow((x-x0),2)+pow((y-y0),2));
    num2 = (R - r)*1;
    return num2;

}

Number function_testreinit (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
  return 0.8*function_two_bubbles (p, parameters,"","");

}

Number function_one_phase (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
    Real num2= -5;
    return num2;
}
Number function_eliptic_bubble2 (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
    Real x = p(0);
    Real y = p(1);
    Real num2=0;
    Real Rx = 0.5;
    Real Ry = 0.525;
    Real x0 = 0;
    Real y0 = 0;
    Real alpha = sqrt(pow((x-x0)/Rx,2)+pow((y-y0)/Ry,2));
    Real theta;
    if (x==0 && y==0)
        theta = -asin(1);
    else
        theta = atan(x/y);

    num2 = (1-alpha)*sqrt(pow(Rx*cos(theta),2)+pow(Ry*sin(theta),2));
    return num2;
}

Number function_supg (const Point& p,
                    const Parameters& parameters,
                    const std::string&,
                    const std::string&)
{
    Real x = p(0);
    Real y = p(1);
    Real x0 = 0;
    Real y0 = 1.5;
    Real R=0.5;
    Real value;

  if ( ((x-x0)*(x-x0)+(y-y0)*(y-y0))<=(R*R))
     value=cos(3.14*(sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/(2*R)));
  else
     value=0;

  return value;
}
Number exact_0_solution(const Point& p,
                      const Parameters&,   // EquationSystem parameters, not needed
                      const std::string&,  // sys_name, not needed
                      const std::string&) // unk_name, not needed);
{
    Real var=0;
    return var;
}

Gradient exact_0_derivative(const Point& p,
                          const Parameters&,   // EquationSystems parameters, not needed
                          const std::string&,  // sys_name, not needed
                          const std::string&) // unk_name, not needed);
{
    Gradient var;
    var.zero();
    return var;
}

RealVectorValue dirichlet_boundary0(const Real x, const Real y, const Real z)
{
    RealVectorValue vector(0,0);
    return vector;
}

RealVectorValue dirichlet_boundary_v_stokes1(const Real x, const Real y, const Real z)
{
    Real U = 1;
    Real u =0;
    Real v = 0;

    u =  -6*U*( Utility::pow<2>(y-0.5)-0.25);

        
    RealVectorValue vector(u,v);

    return vector;
}

Real dirichlet_boundary_p_stokes1(const Real x, const Real y, const Real z){
    Real p=0;
    return p;
}

Real dirichlet_boundary2(const Real x,
                      const Real y,
                      const Real z){

}

bool is_dirichlet_boundary_all(const Real x,
                      const Real y,
                      const Real z){
    return true;
}

bool is_dirichlet_boundary_1point(const Real x,
                      const Real y,
                      const Real z){
    if(x==0 && y==0 && z==0)
        return true;
    else
        return false;
}


bool is_dirichlet_boundary_v_stokes1(const Real x,
                      const Real y,
                      const Real z){
    if((x>7.99999 && x<8.00001) && (y>0.0001 && y < 0.9999))
        return false;
    else
        return true;

    
    
}

bool is_dirichlet_boundary_p_stokes1(const Real x,
                      const Real y,
                      const Real z){
//
    bool flag = is_dirichlet_boundary_v_stokes1(x,y,z);
    if (flag == true)
        return false;
    else
        return true;


}

bool is_dirichlet_boundary_v_stokes2(const Real x,
                      const Real y,
                      const Real z){
    return true;


}

RealVectorValue dirichlet_boundary_v_stokes2(const Real x, const Real y, const Real z)
{
    Real U = 1;
    Real u =0;
    Real v = 0;

    if (y>0.9999)
    {
        u=1;
        if(x>0.95)
        {
            const Real mypi = libMesh::pi;
            u = U - U*sin(mypi/2*(x-0.95)/0.05);
        }
        if(x<0.05)
        {
            const Real mypi = libMesh::pi;
            u = U - U*sin(mypi/2*(0.05-x)/0.05);
        }
    }

    
    RealVectorValue vector(u,v);

    return vector;
}
