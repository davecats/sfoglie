function [sol,prf,flo,tri,blo,eng]=airfoil(prf,flo,tri,blo,eng)

%% PANELIZATION
%  ------------

% Define kinematic viscosity
flo.nu= prf.c*flo.Uinfty/flo.Re;

% Number of wake nodes
eng.NW= round(prf.M/4)+2;   

% x- and y- compoent of incoming flow
flo.ui = flo.Uinfty*cos(flo.alfa); 
flo.vi = flo.Uinfty*sin(flo.alfa);

% Calculate NACA profile nodes
prf = naca4(prf);

% create the panels
prf = create_panels(prf);

%% INVISCID SOLUTION
%  -----------------

% Solve potential flow
[flo]=potential(prf,flo);
% Get stagnation point position
[ prf.Nle,prf.sLE,prf.LE1,prf.LE2 ] = getStagnationPoint( flo.gamma, prf.s );
Nle=prf.Nle;   N=prf.N;
% arc length vector pressure side
prf.sL=prf.s(Nle:end)-(prf.s(Nle)-prf.LE2)*ones(1,N-Nle+1); 
% arc length vector suction side
prf.sU=(prf.s(Nle-1)+prf.LE1)*ones(1,Nle-1)-prf.s(1:Nle-1);   
prf.xU=prf.panels.X(1,1:Nle-1);   prf.xL=prf.panels.X(1,Nle:end);
if prf.sharpTE; prf.xL=[prf.xL,prf.panels.X(2,end) ]; end

%% Create wake
%  --------------

% calculate wake node position , done by integrating ...
%    ... the streamline throug the TE of inviscid solution
flo.wake=GetWakeStreamline(flo,prf,eng.NW);

%% Source influence matrix for airfoil nodes (Bges)
%  ------------------------------------------------

% source influence of airfoil nodes -> i=1,..,N ; j=1,..,N
B=Qlin(prf.nodes.X', prf.nodes.Y' ,prf); % piecewise linear ansatz
% source influence of wake nodes -> i=1,..,N ; j=N+1,..,N+NW
Bw=Qlin(prf.nodes.X', prf.nodes.Y',flo.wake,true); 
Bges=[B, Bw]; Ai=inv(flo.Ages);
% invert airfoil node Coeffs to get Coefficient for Boundary edge velocity U
Btilde=-Ai(1:prf.N,1:prf.N)*Bges;
if prf.sharpTE
    % delete influence of dummy TE panel for sharp TEs
    Btilde(prf.N,prf.N+1:prf.N+flo.wake.N)=0;
end

%% Source influence matrix for wake nodes (C)
%  ------------------------------------------

% influence of airfoil nodes -> i=N+1,..,N+NW ; j=1,..,N+NW
[Cg, Cq] = GradPsiN( prf,flo.wake );
% add gamma influence on Cq
Cq2= Cq - Cg*Btilde;

%% Total source coefficient matrix for EQ Ui= Uinv + Dij qj        XXX
%  --------------------------------------------------------
D= [Btilde;Cq2];
% make sure first wake point has same velocity as trailing edge
D(N+1,:)=D(N,:);
%global arclength vector
sges= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(flo.wake.s)) + flo.wake.s]; 

%% Calculate inviscid velocity on airfoil and wake nodes
%  -----------------------------------------------------

% Velocity at airfoil nodes
UFoil = abs(flo.gamma);
nix= flo.wake.nn(1,:)';%[wake.n(1,1); wake.n(1,:)'];
niy= flo.wake.nn(2,:)';%[wake.n(2,1); wake.n(2,:)'];
UWake = flo.ui*niy - flo.vi*nix + Cg*flo.gamma;
UWake(1)=UFoil(end); % First wake point has same velocity as TE
Uinv=[UFoil; UWake];

%% VISCOUS SOLUTION
%  ----------------

% blowing velocity vector
blo.Vb=zeros(size(Uinv));
% get nodes with blowing and set the blowing velocity vector
blo=addBlowing(blo,prf);
% initial solution guess
sol = GetInitialSolution( prf,flo, tri, eng, Uinv, blo.Vb, 2);
%  coupled boundary layer and potential flow solution
[sol, prf]=NewtonEq( prf,flo,blo,eng,sol,D,Uinv,eng.it);
% If no convergence -> try with different approach for Transition panel EQ
if sol.residual>eng.tol
    eng.tranEQ=true;
    disp('Not converged, Try different Transition approach ')
    disp('------------------------------------------------ ')
    [sol, prf]=NewtonEq( prf,flo,blo,eng,sol,D,Uinv, eng.it );
end

% %% Plots 
% 
% figure 
% hold on
% plot(prf.nodes.X(1:prf.N), sol.Cp(1:prf.N) ,'g')
% title('Pressure coefficient')
% ylabel(' C_p ') 
% xlabel(' x ')
% 
% 
% figure 
% hold on
% plot(prf.nodes.X(1:prf.N), sol.Cf(1:prf.N) ,'g')
% title('Friction coefficient')
% ylabel(' C_f ') 
% xlabel(' x ')
