function [ T, m ,U,DI] = NewtonEq( prf,wake,Uinv,minit,Tinit,D,it)
%NEWTONEQ       sets up the global Newton equationsystem J dz=- f(z) and
%               solves it
%               minit:  initial mass defect vector 
%               Tinit:  initial momentum thickness vector
%               D:      Coefficients of mass Defekt -> U = Uinv + Dm
%               rel:    relaxation factor rel<1 under relaxation, 1<rel<2 over relaxation
%               it:     max number of iteration steps
 

nu=evalin('base','nu');
Nle=prf.Nle;


N=length(Tinit);
T=Tinit;
m=minit;
U=Uinv+D*m;

DI=m./U; %displacement thickness

h=[transpose(prf.panels.L) ;transpose(wake.L)];

% Conditions for closing the system

%Leading edge treatment -> initial condition
% T(s1)= 0.29004*sqrt(s1*nu/U)
T(Nle-1)= 0.29004*sqrt(prf.LE1*nu/U(Nle-1));
T(Nle)  = 0.29004*sqrt(prf.LE2*nu/U(Nle));

% DI(s1)=2.0754*T(s1) -> m(s1)=DI(s1)*U
m(Nle-1)=0.601955021490342*(prf.LE1*nu*U(Nle-1))^0.5;
m(Nle)  =0.601955021490342*(prf.LE2*nu*U(Nle))^0.5;

% first wake node
% Wake: set Quantities for first wake node by means of eq (28)
% Twake_1=Ttop + Tbottom   
T(prf.N+1)=T(1)+T(prf.N);
%     Dwake_1 = Dtop + Dbottom + hTE 
% ->  Dwake_1*Uw1 - Dtop*Uw1 - Dbottom*Uw1 = hTE*Uw1
% ->  mw1 - mtop/Utop*Uw1 - mbottom/Ubottom*Uw1 = hTE*Uw1
m(prf.N+1)=m(1)*U(prf.N+1)/U(1)+m(prf.N)*U(prf.N+1)/U(prf.N);



k=0;
res=1;
while res>1e-5 && k<it
[J,rhs]=JacobiM(T,m,U,DI,D ,h);
dz=J\rhs; % Solution of equation -> gives correction for current Newton-step

[T,m,U,DI]=Update(T,m,Uinv,U,D, dz,prf);


res=max(abs( dz./[T;m]));
k=k+1;
end



% % Set starting conditions to close LGS
% %---------------------------------------------------
% % set m and T to zero at leading edge
% % -> weighted with linear approximation of leading edge position
% % Tin=zeros(1,2*N);
% % Tin(Nle-1)=U(Nle)/(U(Nle)+U(Nle-1));%1/2;
% % Tin(Nle)=U(Nle-1)/(U(Nle)+U(Nle-1));%1/2;
% % min=zeros(1,2*N);
% % min(N+Nle-1)=U(Nle)/(U(Nle)+U(Nle-1));%1/2;
% % min(N+Nle)=U(Nle-1)/(U(Nle)+U(Nle-1));%1/2;
% % 
% % J=[J;Tin;min];
% % rhs=[rhs;0;0];


end

