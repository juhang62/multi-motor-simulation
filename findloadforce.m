function fxs=findloadforce(xs,F,Lstep,ksp)
%%%%%%%%%%%find forces exerted on the motors

k_speff=ksp*Lstep;
m=length(xs);
if m==1
    fxs=F;
else

A=[[-ones(m-1,1) eye(m-1)]; ones(1,m)];
rhs=[xs(2:end)';F/k_speff];
fxs=(A\rhs)*k_speff; %+ hindering; - assisting
end

end