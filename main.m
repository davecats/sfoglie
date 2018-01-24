close all
clear all
%format long eng

addpath('./panel/')
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

Invisc=false;% true;%% only inviscid solution or with Boundary Layer
profile.alfa = 2*pi/180;
profile.Uinfty=1; % Anstroemgeschwindigkeit
Re=6.536*10^4;
%cinematic viscosity
nu= profile.c*profile.Uinfty /Re;%1.53e-5; 

NACA = [4 4 1 2]; % naca profil
NoSkew=true; % if true profile skewness neclected

NW= round(profile.M/4)+2; %   number of wake nodes
%calculates x- and y-component of homogeneous flow
ui = profile.Uinfty*cos(profile.alfa); 
vi = profile.Uinfty*sin(profile.alfa);

% calculate NACA profile Nodes
profile = naca4(profile,NACA,NoSkew);
clear NACA NoSkew

% import profile Nodes
%data=load('Nodes.txt');
%profile.nodes.X=transpose(data(:,1));profile.nodes.Y=transpose(data(:,2));
%profile.N=length(data(:,1)); clear data

% create the panels, identify if Profile has sharp or blunt trailing edge
profile = create_panels(profile);



%%
%  inviscid solution
%------------------------------------

% Solve potential flow
 [field]=potential(profile);  
 
 % get Leading edge position -> secant approximation
 Nle=find(field.gamma<0,1); % first node of pressure side
 profile.Nle=Nle;
% distance from last suction side point
profile.LE1=profile.panels.L(Nle-1)*field.gamma(Nle-1)/(abs(field.gamma(Nle-1))+abs(field.gamma(Nle))); 
% distance from first pressure side point
profile.LE2=profile.panels.L(Nle-1)-profile.LE1;

N=profile.N;
% arc length vector pressure side
su=profile.s(Nle:end)-(profile.s(Nle)-profile.LE2)*ones(1,N-Nle+1); 
% arc length vector suction side
so=(profile.s(Nle-1)+profile.LE1)*ones(1,Nle-1)-profile.s(1:Nle-1); 
profile.so=so;profile.su=su;   

xo=profile.panels.X(1,1:Nle-1);
xu=profile.panels.X(1,Nle:end);
if profile.IsSharp; xu=[xu,profile.panels.X(2,end) ]; end

% Plot polar curves
if Invisc
    Cp=1-field.gamma.^2;
     %sum(Cp(Nle+1:end).*profile.panels.L(Nle:end)')-sum(Cp(1:Nle-1).*profile.panels.L(1:Nle-1)'); %total lift
    if profile.IsSharp; Cp=[Cp; Cp(1)]; end
    figure()
    hold on; box on
    plot([xo,xu(1)], Cp(1:Nle));
    plot(xu,  Cp(Nle:end));
    load naca4412_alfa2.txt
    plot(naca4412_alfa2(:,1),naca4412_alfa2(:,3),'black .');
    xlim([0 1]);
    legend('$C_p=1-\gamma(s)^2$ Saugseite','$C_p=1-\gamma(s)^2$ Druckseite',...
    'Xfoil 160 panels, $\alpha=2^\circ$');
end




%%
%  Wake influence
%------------------------------------



%   calculate wake node position
%------------------------------------
% done by integrating the streamline throug the TE of inviscid solution
wake=GetWakeStreamline(field,profile,NW,1.154);



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
B=field.B; % constant ansatz
%B=Qlin(profile); % piecewise linear ansatz

% influence of wake nodes -> i=1,..,N ; j=N+1,..,N+NW
Bw=Qlin(profile,wake);


Bges=[B, Bw];
Ai=inv(field.A); % invert circulation coefficients
%invert airfoil node Coeffs -> Coefficients by means of eq (10)
Btilde=-Ai*Bges;
% correct the sign for pressure side -> Ui=  gamma_i suction side
%                                    -> Ui= -gamma_i pressure side
%Btilde(Nle:N,:)=-Btilde(Nle:N,:);



% Source coefficient matrix for wake nodes Bges
%-------------------------------------------------------------

% influence of airfoil nodes -> i=N+1,..,N+NW ; j=1,..,N
[Cg, CqFoil] = GradPsi( profile,wake );
Cqw = GradPsiW( wake, 'constant');

 % combine Coefficient matrix for wake nodes by means of eq (31)
Cq=[CqFoil, Cqw];

Cq= Cq + Cg*Btilde;


SourceCoeff= [Bges;Cq];

%  mass defect coefficients 
%---------------------------------------------

sges=[profile.s, wake.s+(wake.s(end)+profile.panels.L(end)/2)*ones(1,NW)]; 
D=MassDefCoeff(SourceCoeff,sges); % forward differences
%D=MassDefCoeff(SourceCoeff,sges,'central'); % central differences


% correct the sign for suction side 
% ->arc length coordinate goes in opposit direction to the stream there
%   q=-dm_ds
D(:,1:Nle-1)=-D(:,1:Nle-1);

% make sure first wake point has same velocity as trailing edge
%D=SourceCoeff;
D(N+1,:)=D(N,:);


%debug
%DT=load('./XFoilWerte/D.txt');
%DT(1:Nle-1,1:Nle-1)=-DT(1:Nle-1,1:Nle-1);
%DT(Nle:end,Nle:end)=-DT(Nle:end,Nle:end);



% Calculate inviscous velocity
%---------------------------------------------

% Velocity at airfoil nodes
%UinvFoil=  Ai*( field.psi0*ones(N,1)+field.t); %-> equal to field.gamma

UFoil = abs(field.gamma);
UWake = ui*wake.nn(2,:).'-vi*wake.nn(1,:).' + Cg*field.gamma;
UWake(1)=UFoil(end); % First wake point has same velocity as TE

CgX=load('./XFoilWerte/Cg.txt');
Uinv=[UFoil; UWake];


clear Do Du Bwake CqFoil Cqwake 
%% 

%  viscous solution
%----------------------------------------------


 
it=3; % maximum number of iterations

%%
% 

% insert the blowing
Vb=zeros(size(Uinv));


 solInit = GetInitialSolution( profile,wake, Uinv,Vb,Re );

%  m=solInit.U.*solInit.D;
%  Unew= Uinv + abs(DT*m);
%  
%  for k=1:it
%     sol = GetInitialSolution( profile,wake, Unew,Vb,Re );
% 
%     m=sol.U.*sol.D;
%     Unew= Uinv + abs(DT*m);
%  end
%  
%  
%  
%  figure
% hold on
% plot(transpose(profile.so),sol.T(1:Nle-1))
% plot(transpose(profile.so),sol.D(1:Nle-1))
% line([profile.so(sol.TranU) profile.so(sol.TranU)], [0 sol.D(sol.TranU)],'color','black');
% title('Saugseite')
% legend('$\delta_2$','$\delta_1$','Transitionspunkt','location','best'); 
%  
% figure
% hold on
% plot(transpose(profile.su),sol.T(Nle:N))
% plot(transpose(profile.su),sol.D(Nle:N))
% line([profile.su(sol.TranL-Nle+1) profile.su(sol.TranL-Nle+1)], [0 sol.D(sol.TranL)],'color','black');
% title('Druckseite')
% legend('$\delta_2$','$\delta_1$','Transitionspunkt','location','best'); 
% 
% figure
% hold on
% plot(transpose(wake.s),sol.T(N+1:N+wake.N))
% plot(transpose(wake.s),sol.D(N+1:N+wake.N))
% title('Nachlauf')
% legend('$\delta_2$','$\delta_1$','$\delta_2$ XFoil','$\delta_1$ XFoil','location','best'); 
 
 
 
 %{
% Newton iterations 
%[T,m,U,d1]=NewtonEq(profile,wake,Uinv,minit,Tinit,DT,it);


% Plots------------------------------------
% pressure side
figure
hold on
plot(xu,d1(Nle:N));
plot(xu,T(Nle:N));
% plot([0;xu],d1(Nle:N+1));
% plot([0;xu],T(Nle:N+1));
plot(xu,Dinit(Nle:N));
plot(xu,Tinit(Nle:N));
title('Druckseite')
legend('$\delta_1$','$\delta_2$','$\delta_1$ Blasius Lsg','$\delta_2$ Blasius Lsg','location','best');


% suction side
%soi=s(Nle)*ones(1,Nle)-s(Nle:-1:1);  
figure
hold on
plot(xo,d1(1:Nle-1));
plot(xo,T(1:Nle-1));
% plot([xo;0],d1(1:Nle));
% plot([xo;0],T(1:Nle));
plot(xo,Dinit(1:Nle-1));
plot(xo,Tinit(1:Nle-1));
title('Saugseite')
legend('$\delta_1$','$\delta_2$','$\delta_1$ Blasius Lsg','$\delta_2$ Blasius Lsg','location','best');


figure
hold on
 plot(wake.x,d1(N+1:end));
 plot(wake.x,T(N+1:end));
% plot(wake.x,d1(N+2:end));
% plot(wake.x,T(N+2:end));
plot(wake.x,Dinit(N+1:end));
plot(wake.x,Tinit(N+1:end));
title('Nachlauf')
legend('$\delta_1$','$\delta_2$','$\delta_1$ Anfangsl\"osung','$\delta_2$ Anfangsl\"osung','location','best');
%}


% "thicked" shape
% d1=Dinit; %tst
% 
% XT=profile.nodes.X - profile.nodes.n(1,:).*transpose(d1(1:profile.N));
% YT=profile.nodes.Y - profile.nodes.n(2,:).*transpose(d1(1:profile.N));
% XTw=wake.x' + wake.nn(1,:).*transpose(d1(profile.N+1:end));
% YTw=wake.y' + wake.nn(2,:).*transpose(d1(profile.N+1:end))/2;
% YTw2=wake.y' - wake.nn(2,:).*transpose(d1(profile.N+1:end))/2;
% 
% figure(); 
% hold on; box on;
% plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
% plot(wake.x,wake.y,'k','Linewidth',1.2);
% plot(XT,YT,'b');
% plot(XTw,YTw,'b');
% plot(XTw,YTw2,'b');
% axis equal; xlabel('x'); ylabel('y') 

%}





