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
%  air foil geometry and panels
%------------------------------------

% Parameters
profile.c = 1;       % scale factor of the profile   
profile.M = 80;      % number of x-values, where nodes will be -> 2*M-1 nodes in total

Invisc=false;% true; % % only inviscid solution or with Boundary Layer
profile.alfa = 0; %5*pi/180;%  

profile.Uinfty=1; % Anstroemgeschwindigkeit
Re= 4e5;  % 6.536*10^4;%
nkrit=  0.15;%  % critical amplification exponent for transition to take place

%cinematic viscosity
nu= profile.c*profile.Uinfty /Re;

NACA = [4 4 1 2]; % naca profil
NoSkew=true; % if true profile skewness neclected

% blowing region
withBlowing=[true;...  % blowing on suction side  
             false];    % blowing on pressure side 
% start region         
xBstart= [0.25;...
          0.25]* profile.c;
% end region      
xBend  = [0.5;...
          0.5]* profile.c;
% blowing intensity      
intensity=[0.001;...
           0.001]* profile.Uinfty;

% tripping
trip=[ false;...       % tripping on suction side 
       false];         % tripping on pressure side 
% trip=[ true;... % tripping on suction side      
%        true];   % tripping on pressure side  
xtrip=[ 0.14;...
        0.08]*profile.c;

%--------------------------------------------------


NW= round(profile.M/4)+2; %   number of wake nodes
%calculates x- and y-component of homogeneous flow
ui = profile.Uinfty*cos(profile.alfa); 
vi = profile.Uinfty*sin(profile.alfa);

% calculate NACA profile Nodes
profile = naca4(profile,NACA,NoSkew);
 clear NACA NoSkew

% % import profile Nodes
% data=load('Nodes.txt');
% profile.nodes.X=transpose(data(:,1));profile.nodes.Y=transpose(data(:,2));
% profile.N=length(data(:,1)); clear data

% create the panels, identify if Profile has sharp or blunt trailing edge
profile = create_panels(profile);



%%
%  inviscid solution
%------------------------------------

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


% Plot polar curves
if Invisc
    Cp=1-field.gamma.^2;
    CL=getCL(profile,field.gamma);
    if profile.IsSharp; Cp=[Cp; Cp(1)]; end
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
wake=GetWakeStreamline(field,profile,NW);


% %Plot wake
% figure(); 
% hold on; box on;
% plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
% plot(wake.x,wake.y,'b');
% axis equal; xlabel('x'); ylabel('y')

if Invisc
   LEstr = LEstreamline( field,profile,round(3*NW/4) ,1.2);
   figure(); 
   hold on; box on;
   plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
   plot(wake.x,wake.y,'b');
   plot(LEstr.x,LEstr.y,'b');
   axis equal; xlabel('x'); ylabel('y') 
   return;
end


% Source coefficient matrix for airfoil nodes Bges
%-------------------------------------------------------------

% influence of airfoil nodes -> i=1,..,N ; j=1,..,N
B=Qlin(profile.nodes.X', profile.nodes.Y' ,profile); % piecewise linear ansatz

% influence of wake nodes -> i=1,..,N ; j=N+1,..,N+NW
Bw=Qlin(profile.nodes.X', profile.nodes.Y' ,wake,true); 

Bges=[B, Bw];


% [L,U,ind] = lu(field.Ages,'vector');
% ALU= U + L-eye(size(field.Ages));

Ai=inv(field.Ages);

%invert airfoil node Coeffs -> Coefficients by means of eq (10)
Btilde=-Ai(1:profile.N,1:profile.N)*Bges;


% Source coefficient matrix for wake nodes Bges
%-------------------------------------------------------------

% influence of airfoil nodes -> i=N+1,..,N+NW ; j=1,..,N
[Cg, Cq] = GradPsiN( profile,wake );


% add gamma influence on Cq
Cq2= Cq - Cg*Btilde;


D= [Btilde;Cq2];
% make sure first wake point has same velocity as trailing edge
D(N+1,:)=D(N,:);
     

%global arclength vector
sges= [profile.s, (profile.s(end)+profile.panels.L(end)/2)*ones(size(wake.s)) + wake.s]; 


% Calculate inviscous velocity
%---------------------------------------------

% Velocity at airfoil nodes
%UinvFoil=  Ai*( field.psi0*ones(N,1)+field.t); %-> equal to field.gamma

UFoil = abs(field.gamma);
nix= wake.nn(1,:)';%[wake.n(1,1); wake.n(1,:)'];
niy= wake.nn(2,:)';%[wake.n(2,1); wake.n(2,:)'];

UWake = ui*niy - vi*nix + Cg*field.gamma;
UWake(1)=UFoil(end); % First wake point has same velocity as TE
clear nix niy

Uinv=[UFoil; UWake];


clear Do Du Bwake CqFoil Cqwake 
%% 

%  viscous solution
%----------------------------------------------

%   without blowing

% initial solution for global Newton Method
Vb=zeros(size(Uinv));
ini = GetInitialSolution( profile,wake, Uinv,Vb,Re, 2, trip, xtrip );
%PlotStuff(profile,wake,ini, 'delta');
%PlotStuff(profile,wake,ini, 'U');


it=16; % maximum number of iterations

%  coupled boundary layer and potential flow solution
[sol, prfE]=NewtonEq( profile,wake,ini,D,Uinv,it);


% plot
% inds=(prfE.Nle-1:-1:1);             % suction side node indizes
% indp=(prfE.Nle  :prfE.N);           % pressure side node indizes
% indw=(prfE.N    :prfE.N + wake.N);  % wake node indizes

% PlotStuff(profile,wake,sol, 'delta');
% PlotStuff(profile,wake,sol, 'U');
% PlotStuff(profile,wake,sol, 'tau');
% PlotStuff(profile,wake,sol, 'Cp');




%   with blowing
%%

if ~withBlowing(1) && ~withBlowing(2); return; end

Vb=zeros(size(Uinv));
if withBlowing(1)
    indB1= find( xU < xBend(1) & xU > xBstart(1));
    Vb(indB1)=intensity(1);
end
if withBlowing(2)
    indB2= find( xL < xBend(2) & xL > xBstart(2));
    indB2= indB2 + (profile.Nle-1)*ones(size(indB2));
    Vb(indB2)=intensity(2);
end

iniB = GetInitialSolution( profile,wake, Uinv,Vb,Re, 2, trip, xtrip  );
%PlotStuff(profile,wake,iniB, 'delta');


%  coupled boundary layer and potential flow solution
[solB, prfB]=NewtonEq( profile,wake,iniB,D,Uinv,it);

% indBs=(prfB.Nle-1:-1:1);            % suction side node indizes
% indBp=(prfB.Nle :prfB.N);           % pressure side node indizes
% indBw=(prfB.N   :prfB.N + wake.N);  % wake node indizes


% PlotStuff(profile,wake,solB, 'delta');
% PlotStuff(profile,wake,solB, 'U');
% PlotStuff(profile,wake,solB, 'tau');
% PlotStuff(profile,wake,solB, 'Cp');

%PlotStuff(prfB,wake,solB, 'tau', indBs);

%%
writeOut=false;

if writeOut
    
    tmp=Re; i=0;
    while tmp>1
        tmp=tmp/10;
        i=i+1;  
    end
    first=round(tmp*10);
    str1=['Re',num2str(first),'e',num2str(i-1)];
    alf=round(profile.alfa*180/pi);
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
            tmp=100*nkrit;
            stmp=num2str(tmp);
            if strcmp(stmp(2),'0'); stmp=stmp(1); end
            str3=['_N0',stmp];
        end

    end

    str=[str1,str2,str3];

    Xges=[profile.nodes.X'; wake.x];
    Yges=[profile.nodes.Y'; wake.y];

    datBlow  = [Xges,Yges,sges', solB.D,solB.T,solB.U,solB.Cp,solB.tau, solB.Cf ];
    datNoBlow= [Xges,Yges,sges', sol.D,sol.T,sol.U,sol.Cp,sol.tau, sol.Cf ];

     dlmwrite(['./ComparisonData/Blowing_',str], datBlow);
     dlmwrite(['./ComparisonData/NoBlowing_',str], datNoBlow);
end
%--------------------------------------------------------------------


% calculate drag reduction coefficients
%--------------------------------------------------------------------

NU= prfE.Nle-1;% find(profile.nodes.Y<0,1) -1;
indU= NU:-1:1;
indL= NU+1:profile.N;

r=zeros(max(NU,profile.N-NU) ,2);

r(1:NU,1)= 1- solB.Cf(indU)./sol.Cf(indU);
r(1:profile.N-NU,2)= 1- solB.Cf(indL)./sol.Cf(indL);

cfMU =  0.5*[sol.Cf(indU(1))  ; ( sol.Cf(indU(2:end))  + sol.Cf(indU(1:end-1))  )];
cfMBU=  0.5*[solB.Cf(indU(1)) ; ( solB.Cf(indU(2:end)) + solB.Cf(indU(1:end-1)) )];

sU =[0; prfE.sU(end:-1:1)']; dsU = sU(2:end) - sU(1:end-1);
sBU=[0; prfB.sU(end:-1:1)']; dsBU= sBU(2:end)-sBU(1:end-1);

cfML =  0.5* [sol.Cf(indL(1))   ;( sol.Cf(indL(2:end))  + sol.Cf(indL(1:end-1))  )];
cfMBL=  0.5* [solB.Cf(indL(1))  ;( solB.Cf(indL(2:end)) + solB.Cf(indL(1:end-1)) )];

sL =[0; prfE.sL']; dsL = sL(2:end) - sL(1:end-1);
sBL=[0; prfB.sL']; dsBL= sBL(2:end)-sBL(1:end-1);

sol.CfintU= cumsum(cfMU.*dsU);
sol.CfintL= cumsum(cfML.*dsL);

solB.CfintU= cumsum(cfMBU.*dsBU);
solB.CfintL= cumsum(cfMBL.*dsBL);


R=zeros(max(NU,profile.N-NU) ,2);

R(1:NU,1)= 1-solB.CfintU./sol.CfintU;
R(1:profile.N-NU,2)=1- solB.CfintL./sol.CfintL;

%%

% comparison plots with blowing and without blowing
%---------------------------------------------------

figure 
hold on
plot(prfE.nodes.X(indU), sol.CfintU ,'g')
plot(prfE.nodes.X(indU), solB.CfintU ,'b')
line([prfE.nodes.X(indB1(1))   prfE.nodes.X(indB1(1))]  , [min(sol.CfintU) max(sol.CfintU)],'color','black');
line([prfE.nodes.X(indB1(end)) prfE.nodes.X(indB1(end))], [min(sol.CfintU) max(sol.CfintU)],'color','black');
title('integral friction coefficient')
ylabel(' [c_f] ') 
xlabel(' x ')
legend('without blowing','with blowing','blowing region','location','northeast'); 

figure 
hold on
plot(prfE.nodes.X(indU), r(:,1) ,'g')
plot(prfE.nodes.X(indU), R(:,1) ,'b')
line([prfE.nodes.X(indB1(1))   prfE.nodes.X(indB1(1))]  , [0 1],'color','black');
line([prfE.nodes.X(indB1(end)) prfE.nodes.X(indB1(end))], [0 1],'color','black');
title('Drag-reduction coefficient')
ylabel(' r ') 
xlabel(' x ')
legend('local r','global R','blowing region','location','northeast'); 



%--------------------------------------------------------------------
figure 
hold on
plot(prfE.nodes.X(1:prfE.N), sol.Cp(1:prfE.N) ,'g')
plot(prfE.nodes.X(1:prfE.N), solB.Cp(1:prfE.N) ,'b')
line([prfE.nodes.X(indB1(1))   prfE.nodes.X(indB1(1))]  , [-1 1],'color','black');
line([prfE.nodes.X(indB1(end)) prfE.nodes.X(indB1(end))], [-1 1],'color','black');
title('Pressure coefficient')
ylabel(' C_p ') 
xlabel(' x ')
legend('without blowing','with blowing','blowing region','location','northeast'); 

figure 
hold on
plot(prfE.nodes.X(1:prfE.N), sol.tau(1:prfE.N) ,'g')
plot(prfE.nodes.X(1:prfE.N), solB.tau(1:prfE.N) ,'b')
line([prfE.nodes.X(indB1(1))   prfE.nodes.X(indB1(1))]  , [min(solB.tau) max(sol.tau)],'color','black');
line([prfE.nodes.X(indB1(end)) prfE.nodes.X(indB1(end))], [min(solB.tau) max(sol.tau)],'color','black');
title('wall shear stress')
ylabel(' \tau_w ') 
xlabel(' x ')
legend('without blowing','with blowing','blowing region','location','northeast'); 


figure 
hold on
plot(prfE.nodes.X(1:prfE.Nle-1), sol.D(1:prfE.Nle-1) ,'k')
plot(prfE.nodes.X(1:prfE.Nle-1), sol.T(1:prfE.Nle-1) ,'b')
plot(prfE.nodes.X(1:prfE.Nle-1), solB.D(1:prfE.Nle-1) ,'r')
plot(prfE.nodes.X(1:prfE.Nle-1), solB.T(1:prfE.Nle-1) ,'g')
title('BL thickness on suction side')
ylabel(' \delta ') 
xlabel(' x ')
legend('$\delta_1$ without blowing','$\delta_2$ without blowing','$\delta_1$ with blowing','$\delta_2$ with blowing','location','northeast'); 

