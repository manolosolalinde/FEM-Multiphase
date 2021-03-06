%Reference Density
Rho = 1000;
%Reference viscosity
Mu = 1e-3;
%Reference velocity
U = 0.001;
%Reference Length
L = 0.001;
%Gravity
g = 0;
%Surface tension
sigma = 0;

Re = Rho*L*U/Mu;
Fr = U / sqrt(g*L);
We = Rho*U^2*L/sigma;

output_p = 11.94;
real_p = Rho*U^2*output_p;

output_x = 1;
real_x = L*output_x;

output_poverx = 1;
real_poverx = Rho*U^2/L*output_poverx;


expected_dpdx = -12*Mu*U/L^2


%filename = 'out_t_a_04.ex2.0080';
%filename = 'out_t_a_03.ex2.0023';
filename = 'out_t_a_02.ex2.0015';

pressure = getnc(filename,'vals_nod_var10',[-1,1],[1,-1])';
coordy = getnc(filename,'coordy');
coordx = getnc(filename,'coordx');
%time = getnc(filename,'time_whole')';




i=1;
while i<=length(coordy)
    if coordy(i)~=0.5
        coordx(i) = [];
        coordy(i) = [];
        pressure(i) = [];
        i=i-1;
    end
    i=i+1;
end

real_pressure = Rho*U^2*pressure';
real_x = L*coordx;
x = real_x;

cfun = fit(real_x,real_pressure,'pchipinterp');%pchipinterp %splineinterp
eps = L*0.001;
%plot(x,(cfun(x+eps)-cfun(x))/eps);
gradient = (cfun(x+eps)-cfun(x))/eps;
myvar = expected_dpdx*ones(1,length(x));
%plot(x,myvar');
%plot(x,cfun(x))
clear ymatrix;
ymatrix = [];
ymatrix(:,1) = gradient;
ymatrix(:,2) = myvar;
createfigure2(x,ymatrix);
xlim([L*1,L*7]);
ylim([min(myvar)*1.05,min(myvar)*0.95]);

suma=0;
conta=0;
for i=1:length(x)
   if ((x(i)>=L*2) && (x(i)<=L*6))
       suma = suma + gradient(i);
       conta = conta+1;
   end
end
average = suma/conta


    
    
    
    
    
    