function  sol  = GetInitialSolution( prf, flo, tri, eng, Uinv, Vb, InVal)
%GETINITIALSOLUTION integrates the Boundary equations with the inviscid velocity using 
%                   the initial values from Eppler 1980:
%                   delta2=0.29004*sqrt( s_1*nu/U_1)
%                   delta1=2.0754 *delta2

wake=flo.wake;
Re=flo.Re;
nu=flo.nu;


N=prf.N;
Nle=prf.Nle;

dim=size(Uinv);

% parameters for solution
sol.HmaxLam=3.8; % maximum H values for seperation to take place
sol.HmaxTurb=2.5;
sol.itmax=25; % maximum Newton iterations (XFoil 25)
sol.resmax=1e-12; % maximum residuum to stop Newton iterations

% Tripping
sol.Tripping=tri.active;
sol.xT=tri.x;
sol.sT=[1;1];% value irrelevant, not used



sol.Ret=zeros(dim);
sol.HK=zeros(dim);
sol.T=zeros(dim);
sol.D=zeros(dim);
sol.c=zeros(dim); % vector with amplification factor or Ctau

sol.Vb=zeros(dim);
sol.Vb(1:length(Vb))=Vb;
sol.U=Uinv;

%----------------------------- Set initial values -----------------------

% initial values for T and D

if InVal==1 % -> Eppler 1980
    % suction side
    Ts= 0.29004*sqrt(prf.LE1*nu/sol.U(Nle-1)); 
    Ds= 2.0754*Ts;

    sol.T(Nle-1)=Ts;
    sol.D(Nle-1)=Ds;
    sol.HK(Nle-1)=2.0754;
    
    %pressure side
    Ts= 0.29004*sqrt(prf.LE2*nu/sol.U(Nle));
    Ds= 2.0754*Ts;
    sol.T(Nle)=Ts;
    sol.D(Nle)=Ds;
    sol.HK(Nle)=2.0754;
else % -> Thwaite
    % suction side
    Ts=sqrt(0.075*prf.LE1/(sol.U(Nle-1)*Re) ); % 
    Ds=2.2*Ts;
    sol.T(Nle-1)=Ts;
    sol.D(Nle-1)=Ds;
    sol.HK(Nle-1)=2.2;    
    
    Ts=sqrt(0.075*prf.LE2/(sol.U(Nle)*Re) );
    Ds=2.2*Ts;
    sol.T(Nle)=Ts;
    sol.D(Nle)=Ds;
    sol.HK(Nle)=2.2;
end

% solves EQ from stagnation point to first node to improve the starting
% Node values
[ sol.T(prf.Nle-1), sol.D(prf.Nle-1) ] = ImproveStartNodeSol( nu, sol.U(prf.Nle-1),prf.LE1,sol.Vb(prf.Nle-1),Ts,Ds);
[ sol.T(prf.Nle)  , sol.D(prf.Nle)   ] = ImproveStartNodeSol( nu, sol.U(prf.Nle)  ,prf.LE2,sol.Vb(prf.Nle)  ,Ts,Ds);

% initial node values
sol.c(prf.Nle-1)=0;
sol.c(prf.Nle)=0;

sol.HK(prf.Nle-1) =sol.D(prf.Nle-1)/sol.T(prf.Nle-1);
sol.Ret(prf.Nle-1)=sol.T(prf.Nle-1)*sol.U(prf.Nle-1)/nu;

sol.HK(prf.Nle) =sol.D(prf.Nle)/sol.T(prf.Nle);
sol.Ret(prf.Nle)=sol.T(prf.Nle)*sol.U(prf.Nle)/nu;


%     separatly integrate the boundary layers of suction and pressure side as well as wake

%  suction side
%------------------------------------------------------------------------------------

 sol = walkBoundary(prf,wake,sol,flo,eng,1);

 

% pressure side
%------------------------------------------------------------------------------------

 sol = walkBoundary(prf,wake,sol,flo,eng,2);


% wake
%------------------------------------------------------------------------------------

% initial conditions for the wake
sol.T(N+1)=sol.T(1)+sol.T(N);
sol.D(N+1)=sol.D(1) + sol.D(N) + prf.gap; 
sol.HK(N+1)=(sol.D(N+1)-prf.gap)/sol.T(N+1); % TE gap does not go into shape parameter
sol.Ret(N+1)=sol.T(N+1)*sol.U(N+1)/nu;

sol = walkBoundary(prf,wake,sol,flo,eng,3);


% mass defect
sol.m=sol.U.*sol.D;


% calculate Lift coefficient
sol.CL=getCL(prf,sol.U,flo.alfa,flo.Uinfty);


% final values of the modeled quantities Cf, CD and HS
indL=(sol.iTran(1)+1:sol.iTran(2)-1); % laminar node indizes
indT=[1:sol.iTran(1) , sol.iTran(2):prf.N]; % wake node indizes
indW=prf.N+1:prf.N+wake.N; % turbulent node indizes

I=zeros(size(sol.T));
sol.Cf=I;
sol.Cf(indL)= 2*CF2lam(sol.HK(indL), sol.Ret(indL));
sol.Cf(indT)= 2*CF2turb(sol.HK(indT), sol.Ret(indT));

sol.HS=I;
sol.HS(indL)=H32lam(sol.HK(indL));
sol.HS(indT)=H32turb(sol.HK(indT),sol.Ret(indT));
sol.HS(indW)=H32turb(sol.HK(indW),sol.Ret(indW));

sol.Us=0.5*sol.HS.*( -1/3 + 1./(0.75*sol.HK) );

sol.CD=I;
sol.CD(indL)=CD2lam(sol.HK(indL), sol.Ret(indL));
sol.CD(indT)=CD2turb(sol.HK(indT),sol.Ret(indT),sol.HS(indT), sol.Us(indT),sol.c(indT),false);
sol.CD(indW)=CD2turb(sol.HK(indW),sol.Ret(indW),sol.HS(indW), sol.Us(indW),sol.c(indW),true);

sol.CD=0.5*sol.HS.*sol.CD;

sol.Cp=1-sol.U.^2;
sol.tau=0.5*sol.Cf .*sol.U.^2;

end

