function  sol  = GetInitialSolution( prf, wake, Uinv, Vb, Re )
%GETINITIALSOLUTION integrates the Boundary equations using the initial values from Eppler 1980
%                   delta2=0.29004*sqrt( s_1*nu/U_1)
%                   delta1=2.0754 *delta2

N=prf.N;
nu=evalin('base','nu');
Nle=prf.Nle;

dim=size(Uinv);

% prescribe parameter solution vectors
sol.CfM=zeros(dim); sol.Cf=zeros(dim); sol.CDM=zeros(dim); sol.CD=zeros(dim); sol.HS=zeros(dim); sol.Ret=zeros(dim);
sol.UQ=zeros(dim); sol.Us=zeros(dim); sol.CtEQ=zeros(dim); sol.Del=zeros(dim);


% ------- XFoil data -----------------
data =importdata('./XFoilWerte/BLVARS.txt');
dat=[data(1:3:end-2,:), data(2:3:end-1,:),data(3:3:end,:)]; dat=dat(:,1:7);

cX=[dat(prf.Nle-1:-1:1,2); dat(prf.Nle:end,2)];
DX=[dat(prf.Nle-1:-1:1,3); dat(prf.Nle:end,3)];
TX=[dat(prf.Nle-1:-1:1,4); dat(prf.Nle:end,4)];
HX=DX./TX;
UX=[dat(prf.Nle-1:-1:1,5); dat(prf.Nle:end,5)];
CFMX=[dat(prf.Nle-1:-1:1,6); dat(prf.Nle:end,6)];
CDX=[dat(prf.Nle-1:-1:1,7); dat(prf.Nle:end,7)];
%-------------------------------------

%abs(sol.D-DX)./DX;


sol.HmaxLam=3.8; % maximum H values for seperation to take place
sol.HmaxTurb=2.5;
sol.itmax=25; % maximum Newton iterations (XFoil 25)
sol.resmax=1e-5;%1e-12; % maximum residuum to stop Newton iterations
sol.Vb=Vb;
sol.U=Uinv;
sol.HK=zeros(dim);

% initial values for T and D
% -> Eppler 1980
Ts= 0.29004*sqrt(prf.LE1*nu/sol.U(Nle-1)); 
Ds= 2.0754*Ts;


% -> Thwaite
% Ts=sqrt(0.075*prf.LE1/(sol.U(Nle-1)*Re) ); % 
% Ds=2.2*Ts;

sol.T=zeros(dim);
sol.T(Nle-1)=Ts;
sol.D=zeros(dim);
sol.D(Nle-1)=Ds;

% -> Eppler 1980
Ts= 0.29004*sqrt(prf.LE2*nu/sol.U(Nle));
Ds= 2.0754*Ts;

% -> Thwaite
% Ts=sqrt(0.075*prf.LE2/(sol.U(Nle)*Re) );
% Ds=2.2*Ts;
sol.T(Nle)=Ts;
sol.D(Nle)=Ds;


sol.c=zeros(dim); % vector with amplification factor or Ctau


%suction side
%------------------------------------------------------------------------------------

 sol = walkBoundary(prf,wake,sol,1);

 
figure
hold on
plot(transpose(prf.so),sol.T(1:Nle-1))
plot(transpose(prf.so),sol.D(1:Nle-1))
% st=transpose(prf.so); st=[st(1:sol.TranU); st(1)-sol.tran.sU ;st(sol.TranU+1:end)];
% Tt=sol.T(1:Nle-1); Tt=[Tt(1:sol.TranU); sol.tran.TU ;Tt(sol.TranU+1:end)];
% Dt=sol.D(1:Nle-1); Dt=[Dt(1:sol.TranU); sol.tran.DU ;Dt(sol.TranU+1:end)];
% plot(st,Tt)
% plot(st,Dt)
line([prf.so(sol.TranU) prf.so(sol.TranU)], [0 sol.D(sol.TranU)],'color','black');
plot(transpose(prf.so),TX(1:Nle-1))
plot(transpose(prf.so),DX(1:Nle-1))
title('Saugseite')
legend('$\delta_2$','$\delta_1$','Transitionspunkt','$\delta_2$ XFoil','$\delta_1$ XFoil','location','best'); 




% pressure side
%------------------------------------------------------------------------------------

 sol = walkBoundary(prf,wake,sol,2);



figure
hold on
plot(transpose(prf.su),sol.T(Nle:N))
plot(transpose(prf.su),sol.D(Nle:N))
line([prf.su(sol.TranL-Nle+1) prf.su(sol.TranL-Nle+1)], [0 sol.D(sol.TranL)],'color','black');
plot(transpose(prf.su),TX(Nle:N))
plot(transpose(prf.su),DX(Nle:N))
title('Druckseite')
legend('$\delta_2$','$\delta_1$','Transitionspunkt','$\delta_2$ XFoil','$\delta_1$ XFoil','location','best'); 



 
% wake
%------------------------------------------------------------------------------------

 
sol.T(N+1)=sol.T(1)+sol.T(N);
%sol.D(N+1)=sol.D(1)+sol.D(N)+prf.panels.L(end);
sol.D(N+1)=sol.D(1)+sol.D(N)+prf.ScrossT*prf.panels.L(end);


sol = walkBoundary(prf,wake,sol,3);

%sol.D(N+1:N+wake.N)=sol.D(N+1:N+wake.N) + wake.GAP';


figure
hold on
plot(transpose(wake.s),sol.T(N+1:N+wake.N))
plot(transpose(wake.s),sol.D(N+1:N+wake.N))
plot(transpose(wake.s),TX(N+1:N+wake.N))
plot(transpose(wake.s),DX(N+1:N+wake.N))
title('Nachlauf')
legend('$\delta_2$','$\delta_1$','$\delta_2$ XFoil','$\delta_1$ XFoil','location','best'); 

end

