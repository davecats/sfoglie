 % Script for airfoil calculation with plots and a comparison blowing/ no blowing 


close all
 clear all
%format long eng

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

%  prf and panels
%  ------------------
prf.naca = [4 4 1 2];         %  NACA 4-digit prf 
prf.noSkew  = true;          %  if true neglects prf skewness
prf.sharpTE = false;          %  if true modifies NACA prf for shart trailing edge
prf.c = 1;                    %  prf chord length 
prf.M = 90;                   %  number of control points for each surface
prf.pmode = 1;                %  node distribution mode: 1 (more nodes in middle) 2 (more nodes at LE and TE)

%  Flow
%  ----
flo.alfa = 5*pi/180;         %  Angle of attack 
flo.invisc = false;          %  only inviscid solution 
flo.Uinfty=1;                %  velocity of outer flow
flo.Re= 1e5;                 %  Chord Reynoldsnumber
flo.nkrit=  0.145;%          %  critical amplification exponent for transition

%  Tripping
%  --------
tri.active=[ false;...         %  tripping on suction side 
             false];           %  tripping on pressure side 

tri.x = [ 0.145;...             %  tripping location on suctoin side
          0.29 ]*prf.c;        %  tripping location on pressure side


% ----------- blowing -----------------------
withBlowing=[true;...  % blowing on suction side  
             false];   % blowing on pressure side 
% blowing region
% startpoint    
xBstart= [0.25;...
          0.25]* prf.c;
% end point     
xBend  = [0.86;...
          0.5]* prf.c;
      
% blowing intensity      
intensity=[0.001;...
           0.001]* flo.Uinfty;

pressureCor=false;% true; % include correction Term for pressure

%  Newton Solver
%  --------------
eng.it=20;                   % maximum number of Newton step iterations
eng.tranEQ=false;            % false) xfoil transition EW, true) modified transition EQ
eng.tol=5e-4;                % tolerance of Newton method


%--------------------------------------------------
%cinematic viscosity
flo.nu= prf.c*flo.Uinfty /flo.Re;

eng.NW= round(prf.M/4)+2; %   number of wake nodes

% x- and y- compoent of incoming flow
flo.ui = flo.Uinfty*cos(flo.alfa); 
flo.vi = flo.Uinfty*sin(flo.alfa);

% calculate NACA profile Nodes
prf = naca4(prf);


% % import profile Nodes
% data=load('e387.txt');
% prf.nodes.X=transpose(data(:,1));prf.nodes.Y=transpose(data(:,2));
% prf.N=length(data(:,1)); clear data


% create the panels, identify if profile has sharp or blunt trailing edge
prf = create_panels(prf);


%%

%----------------------------------------------------------------------------
%  inviscid solution
%----------------------------------------------------------------------------


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


% Plot inviscid pressure coefficient
if flo.invisc
    Cp=1-flo.gamma.^2;
    CL=getCL(prf,flo.gamma);
    if prf.sharpTE; Cp=[Cp; Cp(1)]; end
    figure()
    hold on; box on
    plot([xU,xL(1)], Cp(1:Nle));
    plot(xL,  Cp(Nle:end));
    xlim([0 1]);
    legend('$C_p=1-\gamma(s)^2$ Saugseite','$C_p=1-\gamma(s)^2$ Druckseite');
end




%%
%  Wake influence
%------------------------------------


%   calculate wake node position
%------------------------------------
% done by integrating the streamline throug the TE of inviscid solution
flo.wake=GetWakeStreamline(flo,prf,eng.NW);


% %Plot wake
% figure(); 
% hold on; box on;
% plot([prf.panels.X]',[prf.panels.Y]','k','Linewidth',2);
% plot(flo.wake.x,flo.wake.y,'b');
% axis equal; xlabel('x'); ylabel('y')

% Plot profile and wake
if flo.invisc
   figure; 
   hold on; box on;
   plot([prf.panels.X]',[prf.panels.Y]','k','Linewidth',2);
   plot(flo.wake.x,flo.wake.y,'b');
   axis equal; xlabel('x'); ylabel('y') 
   return;
end


% Source coefficient matrix for airfoil nodes Bges
%-------------------------------------------------------------

% source influence of airfoil nodes -> i=1,..,N ; j=1,..,N
B=Qlin(prf.nodes.X', prf.nodes.Y' ,prf); % piecewise linear ansatz

% source influence of wake nodes -> i=1,..,N ; j=N+1,..,N+NW
Bw=Qlin(prf.nodes.X', prf.nodes.Y' ,flo.wake,true); 

Bges=[B, Bw];


% [L,U,ind] = lu(flo.Ages,'vector');
% ALU= U + L-eye(size(flo.Ages));

Ai=inv(flo.Ages);

% invert airfoil node Coeffs to get Coefficient for Boundary edge velocity U
Btilde=-Ai(1:prf.N,1:prf.N)*Bges;


if prf.sharpTE
    % delete influence of dummy TE panel for sharp TEs
    Btilde(prf.N,prf.N+1:prf.N+flo.wake.N)=0;
end

% Source coefficient matrix for wake nodes C
%-------------------------------------------------------------

% influence of airfoil nodes -> i=N+1,..,N+NW ; j=1,..,N+NW
[Cg, Cq] = GradPsiN( prf,flo.wake );


% add gamma influence on Cq
Cq2= Cq - Cg*Btilde;

% Total source coefficient matrix for EQ Ui= Uinv + Dij qj
%-------------------------------------------------------------

D= [Btilde;Cq2];
% make sure first wake point has same velocity as trailing edge
D(N+1,:)=D(N,:);
     

%global arclength vector
sges= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(flo.wake.s)) + flo.wake.s]; 


% Calculate inviscous velocity on airfoil and wake nodes
%---------------------------------------------

% Velocity at airfoil nodes
%UinvFoil=  Ai*( flo.psi0*ones(N,1)+flo.t); %-> equal to flo.gamma

UFoil = abs(flo.gamma);
nix= flo.wake.nn(1,:)';%[flo.wake.n(1,1); flo.wake.n(1,:)'];
niy= flo.wake.nn(2,:)';%[flo.wake.n(2,1); flo.wake.n(2,:)'];

UWake = flo.ui*niy - flo.vi*nix + Cg*flo.gamma;
UWake(1)=UFoil(end); % First wake point has same velocity as TE
clear nix niy

Uinv=[UFoil; UWake];


%% 

%----------------------------------------------------------------------------
%  viscous solution
%----------------------------------------------------------------------------


%   without blowing
%----------------------------------------------

% initial solution of Boundary Layer for global Newton Method
Vb=zeros(size(Uinv));

% initial solution guess
ini = GetInitialSolution( prf,flo, tri, eng, Uinv, Vb, 2);
%  coupled boundary layer and potential flow solution
[sol, prfE]=NewtonEq( prf,flo,eng,ini,D,Uinv,eng.it,pressureCor);

% If no convergence -> try with different approach for Transition panel EQ
if sol.residual>eng.tol
    eng.tranEQ=true;
    disp('Not converged, Try different Transition approach ')
    disp('------------------------------------------------ ')
    if sol.residual<1
        [solTST, prfTST]=NewtonEq( prfE,flo,eng,sol,D,Uinv,eng.it,pressureCor);
    else
        [solTST, prfTST]=NewtonEq( prf,flo,eng,ini,D,Uinv,eng.it,pressureCor);
    end
    
    if solTST.residual < sol.residual
       sol=solTST; prfE=prfTST;
    end
    clear solTST  prfTST TranEQ2
end
%--------------------------------------------------------------


% plot

% PlotStuff(prfE,flo.wake,sol, 'tau');
% PlotStuff(prfE,flo.wake,sol, 'Cp');

% PlotProfile(prfE,flo.wake,sol, 2);


%%

%   with blowing
%----------------------------------------------

if ~withBlowing(1) && ~withBlowing(2); return; end

disp('With Blowing')
disp('-----------------------')

Vb=zeros(size(Uinv));

% get nodes with blowing and set the blowing velocity vector
if withBlowing(1)
    indB1= find( prf.xU < xBend(1) & prf.xU > xBstart(1));
    Vb(indB1)=intensity(1);
end
if withBlowing(2)
    indB2= find( prf.xL < xBend(2) & prf.xL > xBstart(2));
    indB2= indB2 + (prf.Nle-1)*ones(size(indB2));
    Vb(indB2)=intensity(2);
end

% initial solution guess
iniB = GetInitialSolution( prf,flo, tri, eng, Uinv, Vb, 2);
%  coupled boundary layer and potential flow solution
[solB, prfB]=NewtonEq( prf,flo,eng,iniB,D,Uinv,eng.it,pressureCor);

% If no convergence -> try with different approach for Transition panel EQ
if sol.residual>eng.tol
    eng.tranEQ=true;
    disp('Not converged, Try different Transition approach ')
    disp('------------------------------------------------ ')
    if sol.residual<1
        [solTST, prfTST]=NewtonEq( prfB,flo,eng,solB,D,Uinv,eng.it,pressureCor);
    else
        [solTST, prfTST]=NewtonEq( prf,flo,eng,iniB,D,Uinv,eng.it,pressureCor);
    end
    
    if solTST.residual < sol.residual
       solB=solTST; prfB=prfTST;
    end
    clear solTST  prfTST TranEQ2
end
%--------------------------------------------------------------

% PlotStuff(prfB,flo.wake,solB, 'tau');
% PlotStuff(prfB,flo.wake,solB, 'Cp');

% PlotProfile(prfB,flo.wake,solB, 3);



% % Plots with comparison blowing/ no blowing
% BlowingComparison(prfE,flo.wake,sol,prfB,solB,1);


% % calculate drag reduction coefficients
% %--------------------------------------------------------------------
% 
% %suction side
% sU=[0,prfE.sU(end:-1:1)];
% sBU=[0,prfB.sU(end:-1:1)];
% 
% % Interpolate basic and blowing solution to same arclength
% B_CFI=spline(sBU',[0;solB.tau(prfB.Nle-1:-1:1)],sU');
% rU=1-[1; B_CFI(2:end)./sol.tau(prfE.Nle-1:-1:1)];
% 
% B_CFintI=spline(sBU',solB.CI_U,sU');
% RU=1-B_CFintI./sol.CI_U;
% RU(1)=0;


%%

% Write out 

writeOut=false;

if writeOut
    
    tmp=flo.Re; i=0;
    while tmp>1
        tmp=tmp/10;
        i=i+1;  
    end
    first=round(tmp*10);
    str1=['Re',num2str(first),'e',num2str(i-1)];
    alf=round(flo.alfa*180/pi);
    str2=['_alfa',num2str(alf)];

    if trip(1)
        tr=round(xtrip(1)*100);
        stmp=num2str(tr);
        if strcmp(stmp(2),'0'); stmp=stmp(1); end
        str3=['_Trip0',stmp];
    else
        if nkrit>1
            str3=['_N',num2str(round(nkrit))];
        else
            tmp=100*flo.nkrit;
            stmp=num2str(tmp);
            if strcmp(stmp(2),'0'); stmp=stmp(1); end
            str3=['_N0',stmp];
        end

    end

    str=[str1,str2,str3];

    Xges=[prf.nodes.X'; flo.wake.x];
    Yges=[prf.nodes.Y'; flo.wake.y];

    datBlow  = [Xges,Yges,sges', solB.D,solB.T,solB.U,solB.Cp,solB.tau, solB.Cf ];
    datNoBlow= [Xges,Yges,sges', sol.D,sol.T,sol.U,sol.Cp,sol.tau, sol.Cf ];

     dlmwrite(['./ComparisonData/Blowing_',str], datBlow);
     dlmwrite(['./ComparisonData/NoBlowing_',str], datNoBlow);
end
%--------------------------------------------------------------------






