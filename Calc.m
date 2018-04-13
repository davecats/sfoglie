 % does the calculation for a single case without plots etc -> run in parameters.m script !


%close all
 
addpath('./panel/')
addpath('./Amplification/')
addpath('./geometry/')
addpath('./utilities/')
addpath('./models/')
addpath('./wake/')
addpath('./BoundaryLayer/')
set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');



%%
%  air foil geometry and flow Parameters
%------------------------------------------------------------------------------

% ------ Profile ---------
NACA = p.NACA; 
NoSkew  =p.NoSkew ;%        % if true profile skewness neclected
sharpTE =p.sharpTE;         % Sharp or blunt trailing edge
profile.c = p.c;          % chord length -> scales the profile  
profile.M = p.M;        % number of x-values, where nodes will be -> 2*M-1 nodes in total

% Angle of attack
profile.alfa = p.alfa*pi/180;

 
profile.Uinfty=p.Uinfty;       % velocity of outer flow
Re= p.Re;                % Chord Reynoldsnumber
nkrit=  p.nkrit;%          % critical amplification exponent for transition to take place


% --------- tripping ----------------
trip=p.trip;
%tripping location
xtrip=p.xtrip;

% ----------- blowing -----------------------
withBlowing= p.withBlowing;
  
xBstart= p.xBstart;  
xBend  = p.xBend;  
   
intensity=p.intensity;

pressureCor=p.pressureCor; 

% max iterations
it=p.it;


%--------------------------------------------------
%cinematic viscosity
nu= profile.c*profile.Uinfty /Re;

NW= round(profile.M/4)+2; %   number of wake nodes
%calculates x- and y-component of homogeneous flow
ui = profile.Uinfty*cos(profile.alfa); 
vi = profile.Uinfty*sin(profile.alfa);

% calculate NACA profile Nodes
profile = naca4(profile,NACA,NoSkew,sharpTE);
 clear NACA NoSkew sharpTE

%   data=load('e387.txt');
%   profile.nodes.X=transpose(data(:,1));profile.nodes.Y=transpose(data(:,2));
%   profile.N=length(data(:,1)); clear data
 


profile = create_panels(profile);


%----------------------------------------------------------------------------
%  inviscid solution
%----------------------------------------------------------------------------

% Solve potential flow
[field]=potential(profile);  
 
 % get Leading edge position -> secant approximation
[ profile.Nle,profile.sLE,profile.LE1,profile.LE2 ] = getStagnationPoint( field.gamma, profile.s );
 
Nle=profile.Nle;
N=profile.N;
% arc length vector pressure side
profile.sL=profile.s(Nle:end)-(profile.s(Nle)-profile.LE2)*ones(1,N-Nle+1); 
% arc length vector suction side
profile.sU=(profile.s(Nle-1)+profile.LE1)*ones(1,Nle-1)-profile.s(1:Nle-1);   

xU=profile.panels.X(1,1:Nle-1);
xL=profile.panels.X(1,Nle:end);
if profile.IsSharp; xL=[xL,profile.panels.X(2,end) ]; end


%%
%  Wake influence
%------------------------------------


%   calculate wake node position
%------------------------------------
% done by integrating the streamline throug the TE of inviscid solution
wake=GetWakeStreamline(field,profile,NW);

% Source coefficient matrix for airfoil nodes Bges
%-------------------------------------------------------------

% source influence of airfoil nodes -> i=1,..,N ; j=1,..,N
B=Qlin(profile.nodes.X', profile.nodes.Y' ,profile); % piecewise linear ansatz

% source influence of wake nodes -> i=1,..,N ; j=N+1,..,N+NW
Bw=Qlin(profile.nodes.X', profile.nodes.Y' ,wake,true); 

Bges=[B, Bw];

Ai=inv(field.Ages);

% invert airfoil node Coeffs to get Coefficient for Boundary edge velocity U
Btilde=-Ai(1:profile.N,1:profile.N)*Bges;


if profile.IsSharp
    % delete influence of dummy TE panel for sharp TEs
    Btilde(profile.N,profile.N+1:profile.N+wake.N)=0;
end


% Source coefficient matrix for wake nodes C
%-------------------------------------------------------------

% influence of airfoil nodes -> i=N+1,..,N+NW ; j=1,..,N+NW
[Cg, Cq] = GradPsiN( profile,wake );


% add gamma influence on Cq
Cq2= Cq - Cg*Btilde;

% Total source coefficient matrix for EQ Ui= Uinv + Dij qj
%-------------------------------------------------------------

D= [Btilde;Cq2];
% make sure first wake point has same velocity as trailing edge
D(N+1,:)=D(N,:);
     

%global arclength vector
sges= [profile.s, (profile.s(end)+profile.panels.L(end)/2)*ones(size(wake.s)) + wake.s]; 


% Calculate inviscous velocity on airfoil and wake nodes
%---------------------------------------------

UFoil = abs(field.gamma);

UWake = ui*wake.nn(2,:)' - vi*wake.nn(1,:)' + Cg*field.gamma;
UWake(1)=UFoil(end); % First wake point has same velocity as TE
clear nix niy

Uinv=[UFoil; UWake];

Vb=zeros(size(Uinv));

% get nodes with blowing and set the blowing velocity vector
if withBlowing(1)
    indB1= find( xU < xBend(1) & xU > xBstart(1));
    Vb(indB1)=intensity(1);
end
if withBlowing(2)
    indB2= find( xL < xBend(2) & xL > xBstart(2));
    indB2= indB2 + (profile.Nle-1)*ones(size(indB2));
    Vb(indB2)=intensity(2);
end

% initial solution of Boundary Layer for global Newton Method
ini = GetInitialSolution( profile,wake, Uinv,Vb,Re, 2, trip, xtrip  );


%  coupled boundary layer and potential flow solution
[sol, prfE]=NewtonEq( profile,wake,ini,D,Uinv,it);


% If no convergence -> try with different approach for Transition panel EQ
if sol.residual>5e-4
    TranEQ2=true;
    disp('Not converged, Try different Transition approach ')
    disp('------------------------------------------------ ')
    if sol.residual<1
        [solTST, prfTST]=NewtonEq( prfE,wake,sol,D,Uinv,it);
    else
        [solTST, prfTST]=NewtonEq( profile,wake,ini,D,Uinv,it);
    end
    
    if solTST.residual < sol.residual
       sol=solTST; prfE=prfTST;
    end
    clear solTST  prfTST TranEQ2
end

% 
% figure
% hold on
% plot(prfE.nodes.X,sol.Cp(1:prfE.N))
% plot(prfE.nodes.X,solTST.Cp(1:prfE.N))
% legend('EQ1','EQ2')
% 
% figure
% hold on
% plot(sges,sol.tau)
% plot(sges,solTST.tau)
% legend('EQ1','EQ2')


















