function [T,m,U,DI] = Update( T,m,Uinv,U,D, dz,prf )
%UPDATE updates the solution from one Newton step and underrelaxes if
%necessary

nu=evalin('base','nu');Nle=evalin('base','Nle');
N=length(T);
dT=dz(1:N);
dm=dz(N+1:end);

Un=Uinv+D*m; % total velocity
dU=Un-U;

Rel=GetRelaxationFactor(dT,T);
tmp=GetRelaxationFactor(Rel.*dm,m);
Rel(tmp<Rel)=tmp(tmp<Rel);
tmp=GetRelaxationFactor(Rel.*dU,U);
Rel(tmp<Rel)=tmp(tmp<Rel);


T=T + Rel.*dT;
m=m + Rel.*dm; 
U=U + Rel.*dU;

% new conditions
T(Nle-1)= 0.29004*sqrt(prf.LE1*nu/U(Nle-1));
T(Nle)  = 0.29004*sqrt(prf.LE2*nu/U(Nle));
m(Nle-1)=0.601955021490342*(prf.LE1*nu*U(Nle-1))^0.5;
m(Nle)  =0.601955021490342*(prf.LE2*nu*U(Nle))^0.5;
T(prf.N+1)=T(1)+T(prf.N);
m(prf.N+1)=m(1)*U(prf.N+1)/U(1)+m(prf.N)*U(prf.N+1)/U(prf.N);

DI= m./U;


end

