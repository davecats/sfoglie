close all
clear all
%format long eng

addpath('./panel/')
addpath('./geometry/')
addpath('./utilities/')
addpath('./models/')
addpath('./wake/')
set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');

%%
%  air foil geometry and panels
%------------------------------------

% Parameters
profile.c = 1;       % scale factor of the profile   
profile.N = 80;%140; % number of nodes on each bottom and top side ->2N-1 panels and nodes
                     % node with index i=N has coordinates (0,0)
                     % node with index i=2*N is the same as i=1
NW=40; %   number of wake nodes
profile.alfa = 2*pi/180;
profile.Uinfty=1; % Anstroemgeschwindigkeit
NACA = [4 4 1 2];

Re=6.536*10^4;
nu= profile.c*profile.Uinfty /Re;%1.53e-5; %cinematic viscosity

%calculates x- and y-component of homogeneous flow
ui = profile.Uinfty*cos(profile.alfa); 
vi = profile.Uinfty*sin(profile.alfa);

% Create profile and panels
profile = naca4(profile,NACA);
profile = create_panels(profile);


% Plot profile
% figure(); hold on; box on;
% plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
% axis equal; xlabel('x'); ylabel('y')

clear NACA
%%
%  inviscid solution
%------------------------------------

% Solve potential flow
 [field]=potential(profile);  
 %Nle=round((profile.M-1)/2); %Leading edge node index
 Nle=find(field.gamma<0);Nle=Nle(1); % first node of pressure side
 
 
% figure()
% hold on; box on
% plot(mean(profile.panels.X(1:Nle,:),2), (1-field.gamma(1:Nle).^2));
% plot(mean(profile.panels.X(Nle:end,:),2), (1-field.gamma(Nle:end).^2));
% load naca4412_alfa2.txt
% plot(naca4412_alfa2(:,1),naca4412_alfa2(:,3),'black .');
% xlim([0 1]);
% legend('$C_p=1-\gamma(s)^2$ Saugseite','$C_p=1-\gamma(s)^2$ Druckseite',...
% 'Xfoil 160 panels');




%%
%  viscous solution
%------------------------------------
N=profile.M;

%vector with panel length
L=profile.panels.L;
% Arc length vectors
s=profile.panels.s;

su=s(Nle:end)-s(Nle)*ones(N-Nle+1,1); % arc length vector pressure side
so=s(Nle-1)*ones(Nle-1,1)-s(1:Nle-1); % arc length vector suction side
   
xo=mean(profile.panels.X(1:Nle-1,:),2);
xu=mean(profile.panels.X(Nle:end,:),2);


 %{
%Uinv=abs(field.gamma(:)); % inviscous part of velocity on airfoil
%Ai=inv(field.A); % invert circulation coefficients
%   Coefficient matrix airfoil only (i=1,..,N; j=1,...,N )
%-----------------------------------------
% BtildeFoil=-Ai*field.B;
% %calculate the coefficients of the mass defect m_i from B_ij q_j 
% DFoil=MassDefCoeff(BtildeFoil,transpose(s)); % q with forward diffs
% 
% % correct the sign of D for suction / pressure side
% DFoil=[DFoil(:,1:Nle-1), -DFoil(:,Nle:end)];
% 
% % airfoil velocity without wake influence 
% %  -> Ui=Uinv_i + sum_1^N B_ij q_j = Uinv_i + sum_1^N D_ij m_j
%}



%%

%   calculate wake node position
%------------------------------------
% done by integrating the streamline throug the TE of inviscid solution
wake=GetWakeStreamline(field,profile,NW,1.13);


% % Plot wake
% figure(); 
% hold on; box on;
% plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
% plot(wake.x,wake.y,'k');
% axis equal; xlabel('x'); ylabel('y')


%  equations for wake nodes -> only airfoil influence (i=N+1,...,N+NW; j=1,...,N)
%----------------------------------------------
[Cg, CqFoil]=GradPsi(profile, wake);



%  Influence of wake sources to integral (i=1,...,N+NW; j=N+1,...,N+NW)
%----------------------------------------------
[ Bwake, Cqwake ] = WakeSourceCoeffs( wake,profile );
% Kutta condition, wake source influence on first wake node neglected
Cqwake(1,:)=zeros(1,NW);
 

 % combine Coefficient matrix for airfoil nodes by means of eq (7)f
B=[field.B, Bwake];
 % combine Coefficient matrix for wake nodes by means of eq (31)
Cq=[CqFoil, Cqwake];


Ai=inv(field.A); % invert circulation coefficients
%invert airfoil node Coeffs -> Coefficients by means of eq (10)
Btilde=-Ai*B;
% correct the sign for pressure side -> Ui=  gamma_i suction side
%                                    -> Ui= -gamma_i pressure side
Btilde(:,Nle:1:N)=-Btilde(:,Nle:1:N);
 
 
% complete matrix for q

Bges= [Btilde; CqFoil, Cqwake];
Bges(N+1:1:end,:)= Bges(N+1:1:end,:)+Cg*Btilde;


% % make sure first wake point has same velocity as trailing edge
Bges(N+1,1:1:N)=Bges(N,1:1:N);



%  mass defect coefficients 
%---------------------------------------------

% complete arc length vector with wake
sges=[s; wake.s+(s(end)+L(end)/2)*ones(NW,1)];   


D=MassDefCoeff(Bges,transpose(sges));

% correct the sign for suction side 
% ->arc length coordinate goes in opposit direction to the stream there
%   q=-dm_ds
D(:,1:1:Nle-1)=-D(:,1:1:Nle-1);




% Calculate inviscous velocity
%---------------------------------------------

% Velocity at airfoil nodes
%UinvFoil=  Ai*( field.psi0*ones(N,1)+field.t); %-> equal to field.gamma

UFoil = abs(field.gamma);
UWake = ui*wake.nn(:,2)-vi*wake.nn(:,1) + Cg*abs(field.gamma);

Uinv=[UFoil; UWake];


clear Do Du Bwake CqFoil Cqwake 
%% 

% solution of the coupled problem with Newton method


 % initial solution -> Blasius solution, plate
 a=@(x) sqrt(x*nu);
 
 % momentum thickness
 T1=0.664*a(xo);%*a(so)
 T2=0.664*a(xu);%*a(su)
 TW=0.664*a(wake.x); % wake
 Tinit=[T1; T2;TW]; 
 
 % displacement thickness
 D1=1.7208*a(xo);%*a(so)
 D2=1.7208*a(xu);%*a(su)
 DW=1.7208*a(wake.x); % wake
 Dinit=[D1;D2;DW];
  clear T1 T2 D1 D2
  
 minit=Uinv.*Dinit; % initial mass defect
 

 Uinit=Uinv+D*minit; % resulting initial velocoty (not necessary) 



rel=0.7;% Relaxation factor rel<1 under relaxation, 1<rel<2 overrelaxation
it=1; % maximum number of Newton iterations

% Newton iterations 
z=NewtonEq(Uinv,minit,Tinit,[L;wake.L],D,it,rel);


% end solution
z=z(:,end);
T=z(1:N+NW);
m=z(N+NW+1:end);

U=Uinv+D*m; % total boundary edge velocity
d1=m./U;    % displacement thickness



% Plots------------------------------------
% pressure side
figure
hold on
plot(xu,d1(Nle:N));
plot(xu,T(Nle:N));
plot(xu,Dinit(Nle:N));
plot(xu,Tinit(Nle:N));
title('Druckseite')
legend('$\delta_1$','$\delta_2$','$\delta_1$ Blasius Lsg','$\delta_2$ Blasius Lsg');


% suction side
%soi=s(Nle)*ones(1,Nle)-s(Nle:-1:1);  
figure
hold on
plot(xo,d1(1:Nle-1));
plot(xo,T(1:Nle-1));
plot(xo,Dinit(1:Nle-1));
plot(xo,Tinit(1:Nle-1));
title('Saugseite')
legend('$\delta_1$','$\delta_2$','$\delta_1$ Blasius Lsg','$\delta_2$ Blasius Lsg');


figure
hold on
plot(wake.x,d1(N+1:end));
plot(wake.x,T(N+1:end));
plot(wake.x,Dinit(N+1:end));
plot(wake.x,Tinit(N+1:end));
title('Nachlauf')
legend('$\delta_1$','$\delta_2$','$\delta_1$ Blasius Lsg','$\delta_2$ Blasius Lsg');
%}







