function [ z ] = NewtonEq( Uinv,minit,Tinit,L,D,it,rel)
%NEWTONEQ       sets up the global Newton equationsystem J dz=- f(z) and
%               solves one Newton-step
%               gamma: circulation distribution
%               D: Coefficients of mass Defekt -> U = Uinv + Dm
%               L: Vector with panel length
%               minit: initial mass defect vector 
%               Tinit: initial momentum thickness vector
%               it: max number of iteration steps
%               rel: relaxation factor rel<1 under relaxation, 1<rel<2 overrelaxation

%               dz: correction for solution vektor z=[T1,..,TN,m1,..,mN]^T
%                   -> 1 column From each iteration

 
if nargin==6; rel=1;end;

n=length(Tinit);
T=Tinit;
m=minit;
%z=[T;m];
k=0;
res=1;
while res>1e-2 && k<it
[J,rhs]=JacobiM(Uinv,m,T,L,D);
dz=J\rhs; % Solution of equation -> gives correction for current Newton-step
  
T=T+rel*dz(1:n);
m=m+rel*dz(n+1:end); %correction from current Newton step
res=max(dz./[T;m]);
k=k+1;
%z(1:2*N,k)=[T;m];
end

z=[T;m];




end

