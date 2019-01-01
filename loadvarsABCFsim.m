clear
global kab kba kbc kcb kcf kfc ksf0 ksb0 kdet0 katt dokbt mtot Lstep ksp
kab=4000;kba=2500;kbc=3000;kcb=2000;kcf=4000;kfc=4000;
ksf0=300;ksb0=50;

kdet0=1;katt=3;
dokbt=0.15; %delta/kbt

mtot=5;

Lstep=8; ksp=0.5;

nruns=2000;
FF=-5:4:15;
nFF=length(FF);
Vsim=zeros(nFF,1); Dsim=zeros(nFF,1);
Xgroup=cell(nFF,1);
for i=1:nFF
[Vsim(i), Dsim(i), Xgroup{i}, tmax]=simsemimarkovstd4(FF(i),nruns,kab,kba,kbc,kcb,kcf,kfc, ksf0, ksb0, kdet0, katt, dokbt, mtot, Lstep, ksp);
end
save simdata.mat