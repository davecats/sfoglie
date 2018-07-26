function [prf,flo,CoeffMatrix,Uinv,sges ] = InviscidSolution(prf,flo,eng)
%INVISCIDSOLUTION executes the inviscid part of the solution process
%                   1. creates the profile geometry and calculate normal vectors etc.  
%                   2. calculates the inviscid flow solution (gamma)
%                   3. calculates the wake streamline with inviscid solution of psi
%                   4. calculates all necessary the coefficient matrices of the BEM
%                   5. calculates the inviscid boundary layer edge velocity Uinv



%% Step 1 : airfoil geometry and discretisation
%-------------------------------------------------------------------------

% Define kinematic viscosity
flo.nu= prf.c*flo.Uinfty/flo.Re;

% Number of wake nodes
eng.NW= round(prf.M/4) + 2;   

% x- and y- compoent of incoming flow
flo.ui = flo.Uinfty*cos(flo.alfa); 
flo.vi = flo.Uinfty*sin(flo.alfa);

% Calculate NACA profile nodes
prf = naca4(prf);

% import profile Nodes in case of a given list
% data=load('e387_N160.txt');
% prf.nodes.X=transpose(data(:,1));prf.nodes.Y=transpose(data(:,2));
% prf.N=length(data(:,1)); clear data

% create the panels
prf = create_panels(prf);


%% Step 2 : inviscid flow solution
%-------------------------------------------------------------------------

% Solve potential flow
[flo]=potential(prf,flo);
% Get stagnation point position
[ prf.Nle,prf.sLE,prf.LE1,prf.LE2 ] = getStagnationPoint( flo.gamma, prf.s );

% arc length vector pressure side
prf.sL=prf.s(prf.Nle:end)-(prf.s(prf.Nle)-prf.LE2)*ones(1,prf.N - prf.Nle+1); 
% arc length vector suction side
prf.sU=( prf.s(prf.Nle-1)+ prf.LE1 )*ones(1,prf.Nle-1)-prf.s(1:prf.Nle-1);  
% x vectors for both sides
prf.xU=prf.panels.X(1,1:prf.Nle-1); 
prf.xL=prf.panels.X(1,prf.Nle:end);
if prf.sharpTE; prf.xL=[prf.xL,prf.panels.X(2,end) ]; end

%% Step 3 : calculate the wake streamline
%-------------------------------------------------------------------------

% calculate wake node position , done by integrating 
%  the streamline through the TE of inviscid solution
flo.wake=GetWakeStreamline(flo,prf,eng.NW);


%% Step 4 : calculate coefficient matrices
%-------------------------------------------------------------------------

% % Source influence matrix for airfoil nodes (Bges)
%  ------------------------------------------------

% source influence of airfoil nodes -> i=1,..,N ; j=1,..,N
CoeffMatrix.B=Qlin(prf.nodes.X', prf.nodes.Y' ,prf); % piecewise linear ansatz
% source influence of wake nodes -> i=1,..,N ; j=N + 1,..,N + NW
CoeffMatrix.Bw=Qlin(prf.nodes.X', prf.nodes.Y',flo.wake,true); 
% total B-matrix
CoeffMatrix.Bges=[CoeffMatrix.B, CoeffMatrix.Bw]; 
% inverted airfoil part
CoeffMatrix.Ai=inv(flo.Ages); 
CoeffMatrix.A=flo.A;
% A-matrix with included Kutta-condition
CoeffMatrix.Ages=flo.Ages;
% invert airfoil node Coeffs to get Coefficient for Boundary edge velocity U
CoeffMatrix.Btilde=-CoeffMatrix.Ai(1:prf.N,1:prf.N)*CoeffMatrix.Bges;

if prf.sharpTE
    % delete influence of dummy TE panel for sharp TEs
    CoeffMatrix.Btilde(prf.N,prf.N+1:prf.N+flo.wake.N)=0;
end

% % Source influence matrix for wake nodes (C)
%  ------------------------------------------

% influence of airfoil nodes -> i=N + 1,..,N + NW ; j=1,..,N + NW
[CoeffMatrix.Cg, CoeffMatrix.Cq] = GradPsiN( prf,flo.wake );
% add gamma influence on Cq
CoeffMatrix.Cq2= CoeffMatrix.Cq - CoeffMatrix.Cg*CoeffMatrix.Btilde;

% % Total source coefficient matrix for EQ Ui= Uinv + Dij qj        XXX
%  --------------------------------------------------------
CoeffMatrix.D= [CoeffMatrix.Btilde;CoeffMatrix.Cq2];
% make sure first wake point has same velocity as trailing edge
CoeffMatrix.D(prf.N+1,:)=CoeffMatrix.D(prf.N,:);
%global arclength vector
sges= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(flo.wake.s)) + flo.wake.s]; 


%% step 5: Calculate inviscid velocity on airfoil and wake nodes
%-------------------------------------------------------------------------

% Velocity at airfoil nodes
UFoil = abs(flo.gamma);

%nix= flo.wake.nn(1,:)';%[flo.wake.n(1,1); flo.wake.n(1,:)'];
%niy= flo.wake.nn(2,:)';%[flo.wake.n(2,1); flo.wake.n(2,:)'];
%UWake = flo.ui*niy - flo.vi*nix + CoeffMatrix.Cg*flo.gamma;

% Velocity on wake nodes
UWake = flo.ui*flo.wake.nn(2,:)' - flo.vi*flo.wake.nn(1,:)' + CoeffMatrix.Cg*flo.gamma;
UWake(1)=UFoil(end); % First wake point has same velocity as TE

Uinv=[UFoil; UWake];




end

