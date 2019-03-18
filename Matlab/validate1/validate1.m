L = 0.04;
U = 0.024295061;
sigma = 0.02361;
mu_ref = 8.9e-4 * 100;
rho_ref = 1000;

R = 0.5*L;
Re = L*U*rho_ref/mu_ref;
We = rho_ref*U^2*L/sigma;

% External fluid density
RhoB = 1;
% Internal fluid density
RhoA = 0.1;
% External fluid viscosity
MuB = 1;
% Internal fluid viscosity std=0.5
MuA = 0.1;

filename = 'out_idate1.ex2.0024';

pressure = getnc(filename,'vals_nod_var10',[-1,1],[1,-1])';
coordy = getnc(filename,'coordy');
coordx = getnc(filename,'coordx');

i=1;
while i<=length(coordx)
    if coordx(i)~=0
        coordx(i) = [];
        coordy(i) = [];
        pressure(i) = [];
        i=i-1;
    end;
    i=i+1;
end


cfun = fit(coordy,pressure','pchipinterp');%pchipinterp %splineinterp
computed_pressure = cfun(0)*rho_ref*U^2
expected_pressure = sigma*(1/R)

relative_error = abs(computed_pressure - expected_pressure)/expected_pressure
    
    
    
    
    
    