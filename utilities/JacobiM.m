function [ J, rhs] = JacobiM( Uinv,m,T,L,D )
%JACOBIM        calculates the Jacobimatrix and right hand side for the Newton system 
%               for solution vector z=[T1,..,TN,m1,..,mN]^T 
%               gamma: circulation distribution
%               D: Coefficients of mass Defekt -> U = Uinv + Dm
%               L: Vector with panel length
%               m: mass defect vector from previous step
%               T: momentum thickness vector from previous step


nu=evalin('base','nu');
Nle=evalin('base','Nle');
NW=evalin('base','NW');

N=length(L)+1; % node number on air foil

h=L;
h(N-NW)=h(N-NW)/2; % first wake point is TE panel midpoint -> take half of TE panel

U=Uinv+D*m; % total velocity
%U=Uinv;

DI=m./U; %displacement thickness

%Filter unrealistic DI
% ind=find(DI(1:Nle-1)<1.02*T(1:Nle-1)); %suction side
% DI(ind)=0.02*T(ind);
% ind=find(DI(Nle:end)<1.00005*T(Nle:end)); % pressure side
% ind=ind+(Nle-1)*ones(length(ind),1);
% DI(ind)=1.00005*T(ind);


% panel midpoint approximations (without TE panel)
TM=(T(2:end)+T(1:end-1))/2;
DM=(DI(2:end)+DI(1:end-1))/2;
%mM=(m(2:end)+m(1:end-1))/2;
UM=(U(2:end)+U(1:end-1))/2;

% trailing edge panel
% TTE=(T(1)+T(end))/2;
% DTE=(DI(1)+DI(end))/2;
% mTE=(m(1)+m(end))/2;


% help quantities

H=DI./T; %shape parameter H1
%HM=(H(2:end)+H(1:end-1));%HM=DM./TM;
Ret= U.*T/nu;
%ReM=(Ret(2:end)+Ret(1:end-1));%UM.*TM/nu;


 
% calculate Cf, CD, H32 and their derivates in respect to H and Ret
[ Cf, dCf_dH, dCf_dRet ] = CF2lam( H, Ret);
[ H32, dH32_dH ]=H32lam(H);
[ CD, dCD_dH, dCD_dRet ]=CD2lam(H,Ret,H32);
dCD_dH32=CD./H32;


% mid point values
CfM=(Cf(2:end)+Cf(1:end-1))/2;
CDM=(CD(2:end)+CD(1:end-1))/2;
H32M=(H32(2:end)+H32(1:end-1))/2;

% Derivates T,D, U
%-------------------------------------------------------------------------

%derivates in respect to T
dRet_dT = U/nu;
dH_dT   =-H./T;
dH32_dT = dH32_dH.*dH_dT;

dCf_dT= dCf_dH.*dH_dT + dCf_dRet.*dRet_dT;
dCD_dT= dCD_dH.*dH_dT + dCD_dRet.*dRet_dT + dCD_dH32.*dH32_dT;


%derivates in respect to D
dH_dD   = 1./T;
dCf_dD  = dCf_dH .*dH_dD;
dH32_dD = dH32_dH.*dH_dD;
dCD_dD  = dCD_dH .*dH_dD + dCD_dH32.*dH32_dD;


%derivates in respect to U
dRet_dU= T/nu;
dCf_dU = dCf_dRet.*dRet_dU;
dCD_dU = dCD_dRet.*dRet_dU;

%-------------------------------------------------------------------------

% finite differences of each panel
dT=T(2:end)-T(1:end-1);
dU=U(2:end)-U(1:end-1);
dH32=H32(2:end)-H32(1:end-1);


% momentum equation -> right hand side of Newton System
f1=dT./h + (2*TM + DM).*dU./(h.*UM) - CfM;


% derivates of momentum equation
%-------------------------------------------------------------------------

% momentum thickness T
    % first Term (T2-T1)/h
    da1_dT1= -1./h;
    da1_dT2=  1./h;

    % second Term (2*TM + DM)*(U2-U1)/(h*UM)
    da2_dT1= dU./(h.*UM);
    da2_dT2= dU./(h.*UM);

    % third Term 
    da3_dT1= dCf_dT(1:end-1)/2; % last node is never "1"
    da3_dT2= dCf_dT(2:end)/2;   % first node is never "2"

    df1_dT1= da1_dT1 + da2_dT1 + da3_dT1;
    df1_dT2= da1_dT2 + da2_dT2 + da3_dT2;

    clear da1_dT1 da2_dT1 da3_dT1 da1_dT2 da2_dT2 da3_dT2
% displacement thickness D
    % first Term -> zero
    % second Term (2*TM + DM)*(U2-U1)/(h*UM)
    da2_dD1= dU./(2*h.*UM);
    da2_dD2= dU./(2*h.*UM);

    % third Term
    da3_dD1= dCf_dD(1:end-1)/2; % last node is never "1"
    da3_dD2= dCf_dD(2:end)/2;   % first node is never "2"

    df1_dD1= da2_dD1 + da3_dD1;
    df1_dD2= da2_dD2 + da3_dD2;

    clear  da2_dD1 da3_dD1 da2_dD2 da3_dD2
% tangential velocity U
    % first Term -> zero
    % second Term (2*TM + DM)*(U2-U1)/(h*UM)
    da2_dU1= -(2*TM + DM)./(2*h.*UM) + (2*TM + DM).*dU./(2*h.*UM.^2) ;
    da2_dU2=  (2*TM + DM)./(2*h.*UM) + (2*TM + DM).*dU./(2*h.*UM.^2) ;

    % third Term
    da3_dU1= dCf_dU(1:end-1)/2; % last node is never "1"
    da3_dU2= dCf_dU(2:end)/2;   % first node is never "2"

    df1_dU1= da2_dU1 + da3_dU1;
    df1_dU2= da2_dU2 + da3_dU2;

    clear  da2_dU1 da3_dU1 da2_dU2 da3_dU2
%-------------------------------------------------------------------------


%shape parameter equation -> right hand side of Newton System
f2= TM.*dH32./h + H32M.*(TM-DM).*dU./(h.*UM) - CDM + H32M.*CfM;

% derivates of momentum equation
%-------------------------------------------------------------------------

% momentum thickness T
    % first Term TM*(HS2-HS1)/h
    da1_dT1= dH32./(2*h) - dH32_dT(1:end-1).*TM./h;
    da1_dT2= dH32./(2*h) + dH32_dT(2:end).*TM./h;
    
    % second Term HSM*(TM-DM)*(U2-U1)/(h*UM)
    da2_dT1= dH32_dT(1:end-1).*(TM-DM).*dU./(2*h.*UM)  + H32M.*dU./(2*h.*UM);
    da2_dT2= dH32_dT(2:end).*(TM-DM).*dU./(2*h.*UM)    + H32M.*dU./(2*h.*UM); 

    % third Term -CDM
    da3_dT1= -dCD_dT(1:end-1)/2;
    da3_dT2= -dCD_dT(2:end)/2;

    % fourth Term HSM*CFM
    da4_dT1= dH32_dT(1:end-1).*CfM  + H32M.*dCf_dT(1:end-1)/2;
    da4_dT2= dH32_dT(2:end).*CfM    + H32M.*dCf_dT(2:end)/2;
    
    df2_dT1= da1_dT1 + da2_dT1 + da3_dT1 + da4_dT1;
    df2_dT2= da1_dT2 + da2_dT2 + da3_dT2 + da4_dT2;
    
    clear da1_dT1  da2_dT1  da3_dT1  da4_dT1 da1_dT2  da2_dT2  da3_dT2  da4_dT2
% displacement thickness D
    % first Term TM*(HS2-HS1)/h
    da1_dD1= -dH32_dD(1:end-1).*TM./h;
    da1_dD2=  dH32_dD(2:end).*TM./h;
    
    % second Term HSM*(TM-DM)*(U2-U1)/(h*UM)
    da2_dD1= -H32M.*dU./(2*h.*UM) + dH32_dD(1:end-1).*(TM-DM).*dU./(2*h.*UM);
    da2_dD2= -H32M.*dU./(2*h.*UM) + dH32_dD(2:end).*(TM-DM).*dU./(2*h.*UM);

    % third Term -CDM
    da3_dD1= -dCD_dD(1:end-1);
    da3_dD2= -dCD_dD(2:end);

    % fourth Term HSM*CFM
    da4_dD1= dH32_dD(1:end-1).*CfM  + H32M.*dCf_dD(1:end-1)/2;
    da4_dD2= dH32_dD(2:end).*CfM    + H32M.*dCf_dD(2:end)/2;
    
    df2_dD1= da1_dD1 + da2_dD1 + da3_dD1 + da4_dD1;
    df2_dD2= da1_dD2 + da2_dD2 + da3_dD2 + da4_dD2;

    clear da1_dD1  da2_dD1  da3_dD1  da4_dD1 da1_dD2  da2_dD2  da3_dD2  da4_dD2
% tangential velocity U
    % first Term -> zero
    % second Term HSM*(TM-DM)*(U2-U1)/(h*UM)
    da2_dU1= -H32M.*(TM-DM)./(h.*UM) - H32M.*(TM-DM).*dU./(2*h.*UM.^2);
    da2_dU2=  H32M.*(TM-DM)./(h.*UM) - H32M.*(TM-DM).*dU./(2*h.*UM.^2);

    % third Term -CDM
    da3_dU1= -dCD_dU(1:end-1);
    da3_dU2= -dCD_dU(2:end);

    % fourth Term HSM*CFM
    da4_dU1=  H32M.*dCf_dU(1:end-1)/2;
    da4_dU2=  H32M.*dCf_dU(2:end)/2;
    
    df2_dU1=  da2_dU1 + da3_dU1 + da4_dU1;
    df2_dU2=  da2_dU2 + da3_dU2 + da4_dU2;

    clear  da2_dU1  da3_dU1  da4_dU1  da2_dU2  da3_dU2  da4_dU2

%-------------------------------------------------------------------------


% switch to m as variable

%derivates of D,U in respect to m
tmp=-m./U.^2;

tmp=tmp*ones(1,N);
dD1_dm=tmp(1:end-1,:).*D(1:end-1,:);  % node "1"
dD2_dm=tmp(2:end,:).*D(2:end,:);      % node "2"
% dDi_dmj=1/Ui delta_ij - mi/ui^2 dUi_dmj
dD1_dm=dD1_dm + [diag(1./U(1:end-1)), zeros(N-1,1)];
dD2_dm=dD2_dm + [zeros(N-1,1),diag(1./U(2:end))];

dU1_dm=[diag(1./U(1:end-1)), zeros(N-1,1)];
dU2_dm=[zeros(N-1,1),diag(1./U(2:end))];


% total derivates of equations in respect to m   
df1_dD1=df1_dD1*ones(1,N);
df1_dD2=df1_dD2*ones(1,N);
df1_dU1=df1_dU1*ones(1,N);
df1_dU2=df1_dU2*ones(1,N);

% equation 1-> momentum Eq
df1_dm=  df1_dD1.*dD1_dm + df1_dD2.*dD2_dm ...
        +df1_dU1.*dU1_dm + df1_dU2.*dU2_dm;

df2_dD1=df2_dD1*ones(1,N);
df2_dD2=df2_dD2*ones(1,N);
df2_dU1=df2_dU1*ones(1,N);
df2_dU2=df2_dU2*ones(1,N);

% equation 2-> shape parameter Eq
df2_dm=  df2_dD1.*dD1_dm + df2_dD2.*dD2_dm ...
        +df2_dU1.*dU1_dm + df2_dU2.*dU2_dm;
    

    
% df1_dm= zeros(N-1,N); df2_dm= zeros(N-1,N);
% for k=1:N-1
%      % equation 1-> momentum
%     df1_dm(k,:)=   df1_dD1(k)*dD1_dm(k,:) + df1_dD2(k)*dD2_dm(k,:) ...
%                  + df1_dU1(k)*dU1_dm(k,:) + df1_dU2(k)*dU2_dm(k,:);
%     % equation 2-> shape parameter
%     df2_dm(k,:)=   df2_dD1(k)*dD1_dm(k,:) + df2_dD2(k)*dD2_dm(k,:) ...
%                  + df2_dU1(k)*dU1_dm(k,:) + df2_dU2(k)*dU2_dm(k,:); 
% end


% create Jacobi matrix for Newton System

% contribution of momentum eq
JT1= [diag(df1_dT1), zeros(N-1,1)]; % last node is never a "1" node
JT2= [zeros(N-1,1), diag(df1_dT2)]; % first node is never a "2" node
JTf1= JT1+JT2;
%Jmf1= df1_dm; % mass defect-> all nodes contribute to "1" and "2"


% contribution of shape parameter eq
JT1= [diag(df2_dT1), zeros(N-1,1)]; % last node is never a "1" node
JT2= [zeros(N-1,1), diag(df2_dT2)]; % first node is never a "2" node
JTf2= JT1+JT2;

%Total T part
JT=zeros(2*(N-1),N);
JT(1:2:end-1,1:1:N)=JTf1;
JT(2:2:end,1:1:N)=JTf2;

Jm=zeros(2*(N-1),N);
Jm(1:2:end-1,:)=df1_dm;
Jm(2:2:end,:)=df2_dm;

J=[JT,Jm];

rhs=zeros(2*(N-1),1);
rhs(1:2:end-1)=f1;
rhs(2:2:end)=f2;

%Set starting conditions to close LGS
  
% set TM and mM of LE panel to zero
  Tin=zeros(1,2*N);
  Tin(Nle-1:Nle)=1/2;
  min=zeros(1,2*N);
  min(end,N+Nle-1:N+Nle)=1/2;
    
  J=[J;Tin;min];
  rhs=[rhs;0;0];  
       


end

