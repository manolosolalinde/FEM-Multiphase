%loaddata;

L = 1;%0.04;
U = 1;%0.02429;
sigma = 1;%0.02361;
mu_ref = 9.15826e-2;%8.9e-4 * 10e2;
rho_ref = 1;%1000
% Internal fluid density
RhoA = 100;
% External fluid density
RhoB = 1;
% Internal fluid viscosity std=0.5
MuA = 0.35;%35;
% External fluid viscosity
MuB = 0.1;

Re = L*U*rho_ref/mu_ref;
We = rho_ref*U^2*L/sigma;

Rx = 0.5;
Ry = 0.525; 
%we are ploting the vertical y axis

R = (Ry+Rx)/2;
R = sqrt(Rx*Ry);
D = R*2;
A0 = Ry-Rx;

R2 = R;
n=2;
omega = sqrt((n^3-n)*sigma/((RhoA+RhoB)*rho_ref*R2^3));
tau = R^2*(RhoA+RhoB)*rho_ref/(2*n*((n-1)*MuA*mu_ref+(n+1)*MuB*mu_ref));
An = A0*exp(-time/tau);


%omega = sqrt(24*sigma/((3*RhoA+2*RhoB)*rho_ref*R^3));

analitic = D+(An).*cos(omega*time);


createfigure(time,Data,analitic);


    
    
    
   