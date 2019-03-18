clear;

filename = 'out_anim_a1.ex2'; %13 %8
filename = 'out_anim_idate3_0105b.ex2' %3 %9

level = getnc(filename,'vals_nod_var3',[-1,1],[1,-1])';
velocity = getnc(filename,'vals_nod_var9',[-1,1],[1,-1])';
coordy = getnc(filename,'coordy');
coordx = getnc(filename,'coordx');
time_a = getnc(filename,'time_whole')';

i = size(level,2);
level(:,i) = [];
velocity(:,i) = [];
time_a(i) = [];


i=1;
while i<=length(coordx)
    if coordx(i)~=0
        coordx(i) = [];
        coordy(i) = [];
        level(i,:) = [];
        velocity(i,:) = [];
        i=i-1;
    end;
    i=i+1;
end


for i=1:(length(level))
    cfun = fit(coordy,-level(:,i),'splineinterp');%pchipinterp %splineinterp
    cfun2 = fit(coordy,velocity(:,i),'pchipinterp');%pchipinterp %splineinterp
    x = fminbnd(cfun,min(coordy),max(coordy));
    X(i) = x;
    Curve_a(i) = cfun2(x);
end;

plot(time_a,Curve_a);

%createfigure(time_a,Curve_a);


    
    
    
    
    
    