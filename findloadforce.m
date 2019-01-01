function fxs=findloadforce(xs,F,Lstep,ksp)
% %%%%%%%%%%parameter
% Lstep=8; ksp=0.5;
% %%%%%%
%global Lstep ksp %slow

k_speff=ksp*Lstep;
m=length(xs);
if m==1
    fxs=F;
else
%%%%% xs=[x1-x1 x2-x1 x3-x2]; other codes not written for this 
% tmp=[-1 zeros(1,m-2) 1];
% tmp=repmat(tmp,1,m-1);
% tmp=reshape(tmp,m-1,m);
% A=[tmp;ones(1,m)];

A=[[-ones(m-1,1) eye(m-1)]; ones(1,m)];
rhs=[xs(2:end)';F/k_speff];
fxs=(A\rhs)*k_speff; %+ hindering; - assisting
end

end