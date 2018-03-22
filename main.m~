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

% ------ Profile ---------
NACA = [4 4 1 2]; % naca profil
NoSkew  =false;         % if true profile skewness neclected
sharpTE =false;         % Sharp or blunt trailing edge
profile.c = 1;          % chord length -> scales the profile  
profile.M = 90;        % number of x-values, where nodes will be -> 2*M-1 nodes in total

% Angle of attack
profile.alfa = 5*pi/180;
Invisc = false;% true; %    % only inviscid solution or with Boundary Layer
 

profile.Uinfty=1;       % velocity of outer flow
Re= 4e5;  % 6.536*10^4; % Chord Reynoldsnumber
nkrit=  9;%             % critical amplification exponent for transition to take place


% --------- tripping ----------------
% trip=[ false;...       % tripping on suction side 
%        false];         % tripping on pressure side 
 trip=[ true;... % tripping on suction side      
        true];   % tripping on pressure side  

%tripping location
xtrip=[ 0.14;...
        0.08 ]*profile.c;


% ----------- blowing -----------------------
withBlowing=[true;...  % blowing on suction side  
             false];    % blowing on pressure side 
% blowing region
% startpoint    
xBstart= [0.25;...
          0.25]* profile.c;
% end point     
xBend  = [0.5;...
          0.5]* profile.c;
      
% blowing intensity      
intensity=[0.001;...
           0.001]* profile.Uinfty;



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

% % import profile Nodes
% data=load('Nodes.txt');
% profile.nodes.X=transpose(data(:,1));profile.nodes.Y=transpose(data(:,2));
% profile.N=length(data(:,1)); clear data


% create the panels, identify if profile has sharp or blunt trailing edge
profile = create_panels(profile);


%%

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


% Plot inviscid pressure coefficient
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

% Plot profile and wake
if Invisc
   figure(); 
   hold on; box on;
   plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
   plot(wake.x,wake.y,'b');
   %plot(LEstr.x,LEstr.y,'b');
   axis equal; xlabel('x'); ylabel('y') 
   return;
end


% Source coefficient matrix for airfoil nodes Bges
%-------------------------------------------------------------

% source influence of airfoil nodes -> i=1,..,N ; j=1,..,N
B=Qlin(profile.nodes.X', profile.nodes.Y' ,profile); % piecewise linear ansatz

% source influence of wake nodes -> i=1,..,N ; j=N+1,..,N+NW
Bw=Qlin(profile.nodes.X', profile.nodes.Y' ,wake,true); 

Bges=[B, Bw];


% [L,U,ind] = lu(field.Ages,'vector');
% ALU= U + L-eye(size(field.Ages));

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

% Velocity at airfoil nodes
%UinvFoil=  Ai*( field.psi0*ones(N,1)+field.t); %-> equal to field.gamma

UFoil = abs(field.gamma);
nix= wake.nn(1,:)';%[wake.n(1,1); wake.n(1,:)'];
niy= wake.nn(2,:)';%[wake.n(2,1); wake.n(2,:)'];

UWake = ui*niy - vi*nix + Cg*field.gamma;
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

% PlotStuff(prfE,wake,sol, 'delta');
% PlotStuff(prfE,wake,sol, 'U');
% PlotStuff(prfE,wake,sol, 'tau');
% PlotStuff(prfE,wake,sol, 'Cp');



%%

%   with blowing
%----------------------------------------------

if ~withBlowing(1) && ~withBlowing(2); return; end

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
iniB = GetInitialSolution( profile,wake, Uinv,Vb,Re, 2, trip, xtrip  );
%PlotStuff(profile,wake,iniB, 'delta');


%  coupled boundary layer and potential flow solution
[solB, prfB]=NewtonEq( profile,wake,iniB,D,Uinv,it);

% indBs=(prfB.Nle-1:-1:1);            % suction side node indizes
% indBp=(prfB.Nle :prfB.N);           % pressure side node indizes
% indBw=(prfB.N   :prfB.N + wake.N);  % wake node indizes


% PlotStuff(prfB,wake,solB, 'delta');
% PlotStuff(prfB,wake,solB, 'U');
% PlotStuff(prfB,wake,solB, 'tau');
% PlotStuff(prfB,wake,solB, 'Cp');

%PlotStuff(prfB,wake,solB, 'tau', indBs);

%%

% Write out 

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

%dS=prfB.sLE-prfE.sLE;

%suction side
sU=[0,prfE.sU(end:-1:1)];
sBU=[0,prfB.sU(end:-1:1)];

% Interpolate basic and blowing solution to same arclength
B_CFI=spline(sBU',[0;solB.tau(prfB.Nle-1:-1:1)],sU');
rU=1-[1; B_CFI(2:end)./sol.tau(prfE.Nle-1:-1:1)];

B_CFintI=spline(sBU',solB.CI_U,sU');
RU=1-B_CFintI./sol.CI_U;
RU(1)=0;

%%

% comparison plots with blowing and without blowing
%---------------------------------------------------
DoPlots= false;

if ~DoPlots; return; end

figure 
hold on
plot(sU, sol.CI_U ,'g')
plot(sU, B_CFintI,'b')
line([prfE.sU(indB1(1))   prfE.sU(indB1(1))]  , [min(sol.CI_U) max(sol.CI_U)],'color','black');
line([prfE.sU(indB1(end)) prfE.sU(indB1(end))], [min(sol.CI_U) max(sol.CI_U)],'color','black');
title('integral friction coefficient')
ylabel(' [c_f] ') 
xlabel(' s ')
legend('without blowing','with blowing','blowing region','location','southeast'); 

Cd_red=  (solB.Cdrag-sol.Cdrag)/sol.Cdrag ;
Cnu_red=  (solB.Cnu-sol.Cnu)/sol.Cnu ;
CL_red=  (solB.CL-sol.CL)/sol.CL ;

str={['\Delta C_d=',num2str(100*Cd_red),'%'],['\Delta C_\nu=',num2str(100*Cnu_red),'%'],['\Delta C_L=',num2str(100*CL_red),'%']};

figure 
hold on
plot(sU, rU ,'g')
plot(sU, RU ,'b')
line([prfE.sU(indB1(1))   prfE.sU(indB1(1))]  , [0  max(rU)+0.03],'color','black');
line([prfE.sU(indB1(end)) prfE.sU(indB1(end))], [0  max(rU)+0.03],'color','black');
text(0.8,0.8*max(rU),str);
title('Drag-reduction coefficient')
ylabel(' r ') 
xlabel(' s ')
ylim([min(rU)  max(rU)+0.03])
legend('local r','global R','blowing region','location','northeast'); 



%--------------------------------------------------------------------
figure 
hold on
plot(prfE.nodes.X(1:prfE.N), sol.Cp(1:prfE.N) ,'g')
plot(prfE.nodes.X(1:prfE.N), solB.Cp(1:prfE.N) ,'b')
line([prfE.nodes.X(indB1(1))   prfE.nodes.X(indB1(1))]  , [min(sol.Cp) max(sol.Cp)],'color','black');
line([prfE.nodes.X(indB1(end)) prfE.nodes.X(indB1(end))], [min(sol.Cp) max(sol.Cp)],'color','black');
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






