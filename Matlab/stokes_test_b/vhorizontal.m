importfile('v-horizontal.txt');
%filename = 'out_t_a_04.ex2.0080';
%filename = 'out_t_a_03.ex2.0023';
filename = 'out_t_b_03.ex2.6749';
filename2 = 'out_anim_t_b_03.ex2';

%step = 900;
step = -1;
u = getnc(filename,'vals_nod_var9',[step,-1],[step,-1])';
coordy = getnc(filename,'coordy');
coordx = getnc(filename,'coordx');
%time = getnc(filename,'time_whole')';


i=1;
while i<=length(coordx)
    if coordy(i)~=0.5
        coordx(i) = [];
        coordy(i) = [];
        u(i) = [];
        i=i-1;
    end
    i=i+1;
end

real_u = u';
real_x = coordx;

cfun = fit(real_x,real_u,'pchipinterp');%pchipinterp %splineinterp

count = 0;
for i=0:0.001:1
    count = count +1;
   yvar(count) = i;
   vvar(count) = cfun(i);
end

%plot(x,x10000);%, 'DisplayName', 'x10000 vs y', 'XDataSource', 'y', 'YDataSource', 'x10000');
%hold all;
%plot(cfun);
plot_horizontalv(yvar,vvar,x,x10000);

