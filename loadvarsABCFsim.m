clear
global kab kba kbc kcb kcf kfc ksf0 ksb0 kdet0 katt dokbt mtot Lstep ksp

%%%%%%%%%kinetic parameters
kab=4000;kba=2500;kbc=3000;kcb=2000;kcf=4000;kfc=4000;

%%%%force dependent rate at zero force
ksf0=300;ksb0=50;

kdet0=1; %rate of detachment at zero force
katt=3; %rate of attachment
dokbt=0.15; %delta/kbt

mtot=5; %number of motors

Lstep=8; %length of a phyiscal step
ksp=0.5; %spring constant

nruns=100; %number of runs
FF=-5:4:15; %loading forces 
nFF=length(FF);
Vsim=zeros(nFF,1); %asymptotic velocity
Dsim=zeros(nFF,1); %asymptotic diffusivity
Xgroup=cell(nFF,1); %position of the group 
for i=1:nFF
[Vsim(i), Dsim(i), Xgroup{i}, tmax]=simsemimarkovstd4(FF(i),nruns,kab,kba,kbc,kcb,kcf,kfc, ksf0, ksb0, kdet0, katt, dokbt, mtot, Lstep, ksp);
end
save simdata.mat