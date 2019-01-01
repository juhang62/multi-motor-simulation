function [Vsim, Dsim, Xgroup, tmax]=simsemimarkovstd4(F,nruns,kab,kba,kbc,kcb,kcf,kfc, ksf0, ksb0, kdet0, katt, dokbt, mtot, Lstep, ksp)

tmax=10;
Xgroup=zeros(nruns,1);
parfor irun=1:nruns
t=0;
chst=ones(mtot,1); %make it vertical
xs=zeros(1,mtot);
%i=0;
A0=[0 kab 0 0; kba 0 kbc 0 ; 0 kcb 0 kcf; 0 0 kfc 0; 0 0 0 0];
while t<tmax
    %i=i+1;

chemrates=A0(chst,:);

logatt=~isnan(xs);
xatt=xs(logatt);
xatt=xatt-xatt(1);
fxs=zeros(mtot,0);
fxs(logatt)=findloadforce(xatt,F,Lstep,ksp);
out=zeros(mtot,4);
logA = (chst==1);
logF = (chst==4);
logAorF = any([logA logF],2);
out(logA,1)=ksb0*exp(-fxs(logA)*dokbt.*double(fxs(logA)>0));
out(logF,2)=ksf0*exp(-fxs(logF)*dokbt.*double(fxs(logF)>0));
if sum(logatt)>1
    belldet=min(exp(abs(fxs(logAorF))*dokbt),10); %have to cap kdet; need be same in semi
    out(logAorF,3)=kdet0*belldet;
end
out(chst==5,4)=katt;
phyrates=out;
chmphyrates=[chemrates phyrates];
rates=sum(chmphyrates,2);
t=t+exprnd(1/sum(rates));
actedmotor=randsample(mtot,1,true,rates); %which motor or attach
if isnan(xs(actedmotor)) %attach
    xs(actedmotor)=mean(xs(~isnan(xs)));
    chst(actedmotor)=1;
else
    rates2=chmphyrates(actedmotor,1:7);
    evt2=randsample(7,1,true,rates2); %which reaction    
    if evt2<=4 %pure chemcial
        chst(actedmotor)=evt2;
    elseif evt2==5 %step back
        chst(actedmotor)=4;
        xs(actedmotor)=xs(actedmotor)-1;
    elseif evt2==6 %step forward
        chst(actedmotor)=1;
        xs(actedmotor)=xs(actedmotor)+1;
    else %dettach
        chst(actedmotor)=5; %!! detachment named 5
        xs(actedmotor)=NaN;
    end
end
end
Xgroup(irun)=mean(xs(~isnan(xs)));
end
%Xgroup
Vsim=Lstep*mean(Xgroup)/tmax;
Dsim=Lstep^2*var(Xgroup)/tmax/2;
end
