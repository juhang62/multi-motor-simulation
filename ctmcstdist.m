function out=ctmcstdist(Q)
n=length(Q);
Q(:,end)=ones(n,1);
b=zeros(n,1);
b(end)=1;
out=Q'\b;
end