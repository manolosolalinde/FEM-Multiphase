filename = 'out_anim.ex2';
filename = 'out_anim_output_Wed_Oct_20_150105_2010.ex2';
filename = 'out_anim_4_03.ex2'; %4
filename = 'out_anim_idate2_ereinii_003.ex2';
filename = 'out_anim_idate2_ereinit_0001.ex2';
filename = 'out_anim_idate2_06_03.ex2';%3
filename = 'out_anim_ate2_05.ex2' %4
filename = 'out_anim_date2b_2604a.ex2' %7
filename = 'out_anim_date2b_2604b.ex2' %3
filename = 'out_anim_date2b_2704a.ex2' %3
filename = 'out_anim_date2b_2504.ex2' %7
filename = 'out_anim_date2b_2604a.ex2' %7
filename = 'out_anim_date2b_2804.ex2' %3
filename = 'out_anim_date2b.ex2' %3


level = getnc(filename,'vals_nod_var3',[-1,1],[1,-1])';
coordy = getnc(filename,'coordy');
coordx = getnc(filename,'coordx');
time = getnc(filename,'time_whole')';
size_level = size(level);
i = size_level(2);
level(:,i) = [];
time(i) = [];

i=1;
while i<=length(coordx)
    if coordx(i)~=0
        coordx(i) = [];
        coordy(i) = [];
        level(i,:) = [];
        i=i-1;
    end;
    i=i+1;
end

clear Data;
for i=1:length(time)
    cfun = fit(coordy,level(:,i),'pchipinterp');%pchipinterp %splineinterp
    Data(i) = fzero(cfun,0.5) - fzero(cfun,-0.5);
end;