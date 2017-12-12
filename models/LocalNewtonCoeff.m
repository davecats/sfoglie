function [ JT, JM, f1, f2 ] = LocalNewtonCoeff(i,h,T,Uinv,m,D)
%LOCALNEWTONCOEFF   calculates Jacobi-matrix and right hand sight of each
%                   panel Newton step -> J1 dz1 + J2 dz2 = -f(z1,z2) .
%                   solutionvektor is the triplet [delta_2, delta_1, U]
%                   z1 -> on start node of panel, z2 -> on end node of panel



nu=evalin('base','nu');
% momentum thickness
T1=T(1);
T2=T(1);
% bouandary layer edge velocity
U1=D(1,:)*m +Uinv(1);
U2=D(2,:)*m +Uinv(2);

m1=m(i);
m2=m(i+1);

% displacement thickness
D1=m1/U1;
D2=m2/U2;



% panel midpoint values approximation
TM=(T1+T2)/2;
UM=(U1+U2)/2;
DM=(D1+D2)/2;

% help variables
H1=D1/T1;
H2=D2/T2;
HM=(H1+H2)/2;
Re1=T1*U1/nu;
Re2=T2*U2/nu;
ReM=(Re1+Re2)/2;

%derivates of help variables in respect to T,U,D
%node 1
dH2_dD2=1/T2;
dH2_dT2=-H2/T2;
dRe2_dT2=U2/nu;
dRe2_dU2=T2/nu;
%node 2
dH1_dD1=1/T1;
dH1_dT1=-H1/T1;
dRe1_dT1=U1/nu;
dRe1_dU1=T1/nu;


%  model for Cf/2
%------------------------------------
%node 1
[Cf1, dCf1_dH, dCf1_dRe]=CF2lam(H1,Re1);

% derivates of primary variables
dCf1_dT1 = dCf1_dH*dH1_dT1 + dCf1_dRe*dRe1_dT1;
dCf1_dD1 = dCf1_dH*dH1_dD1                    ;
dCf1_dU1 =                   dCf1_dRe*dRe1_dU1;

%node 2
[Cf2, dCf2_dH, dCf2_dRe]=CF2lam(H2,Re2);

% derivates of primary variables
dCf2_dT2 = dCf2_dH*dH2_dT2 + dCf2_dRe*dRe2_dT2;
dCf2_dD2 = dCf2_dH*dH2_dD2                    ;
dCf2_dU2 =                   dCf2_dRe*dRe2_dU2;

CfM=(Cf1+Cf2)/2;
%correction
%[CfM2, dCfM_dH, dCfM_dRe]=CF2lam(HM,ReM);
%CfM=(CfM+CfM2)/2



%  model for cinetic energy thickness (E)
%------------------------------------
%node 1
[E1,dE1_dH]=delta3lam(T1,H1);

% derivates of primary variables
dE1_dT1= dE1_dH*dH1_dT1 + E1/T1;
dE1_dD1= dE1_dH*dH1_dD1;

%node 2
[E2,dE2_dH]=delta3lam(T2,H2);

% derivates of primary variables
dE2_dT2= dE2_dH*dH2_dT2 + E2/T2;
dE2_dD2= dE2_dH*dH2_dD2;

EM=(E1+E2)/2;

%  model for dissipation integral 2*CD
%------------------------------------

%node 1
[CD1, dCD1_dH, dCD1_dRe]=CD2lam(H1,Re1,E1/T1);

% derivates of prime variables
dHS1_dT1=-E1/T1 + dE1_dT1; % dHs /T1 -> unten schon verrechnet
dHS1_dD1= dE1_dD1;

dCD1_dT1= dCD1_dH*dH1_dT1 + dCD1_dRe*dRe1_dT1 + dHS1_dT1*CD1/E1;
dCD1_dD1= dCD1_dH*dH1_dD1 +                     dHS1_dD1*CD1/E1;
dCD1_dU1=                   dCD1_dRe*dRe1_dU1                  ;


%node 2
[CD2, dCD2_dH, dCD2_dRe]=CD2lam(H2,Re2,E2/T2);

% derivates of prime variables
dHS2_dT2=-E2/T2 + dE2_dT2; % dHs /T1 -> unten schon verrechnet
dHS2_dD2= dE2_dD2;

dCD2_dT2= dCD2_dH*dH2_dT2 + dCD2_dRe*dRe2_dT2 + dHS2_dT2*CD2/E2;
dCD2_dD2= dCD2_dH*dH2_dD2 +                     dHS2_dD2*CD2/E2;
dCD2_dU2=                   dCD2_dRe*dRe2_dU2                  ;

CDM=(CD1+CD2)/2;

%  momentum Equation
%------------------------------------

f1=(T2-T1)/h+(2*TM+DM)*1/UM*(U2-U1)/h-CfM;

% derivates of primary variables of momentum Eq
df1_dT1=-1/h + (U2-U1)/(UM*h) - dCf1_dT1;
df1_dT2= 1/h + (U2-U1)/(UM*h) - dCf2_dT2;
df1_dD1= (U2-U1)/(2*h*UM) - dCf1_dD1;
df1_dD2= (U2-U1)/(2*h*UM) - dCf2_dD2;
df1_dU1= -(2*TM+DM)/(h*UM) -(2*TM+DM)*(U2-U1)/(2*h*UM^2) - dCf1_dU1;
df1_dU2=  (2*TM+DM)/(h*UM) -(2*TM+DM)*(U2-U1)/(2*h*UM^2) - dCf2_dU2;



%  shape parameter Equation
%------------------------------------

f2= (E2-E1)/h - EM/TM* (T2-T1)/h + EM*(1+DM/TM)*(U2-U1)/(h*UM) -CDM + EM/TM*CfM;

% derivates of primary variables of shape parameter Eq
df2_dT1=    EM/(TM*h) + EM/TM^2*(T2-T1)/(2*h) - EM*HM/TM*(U2-U1)/(2*h*UM) ...
          - dE1_dT1/h + dE1_dT1*(T2-T1)/(2*TM*h) + dE1_dT1*(1+DM/TM)*(U2-U1)/(2*h*UM) ...
          + dE1_dT1/(2*TM)*CfM - dCD1_dT1/2 + EM/(2*TM)*dCf1_dT1 - EM/(2*TM^2)*CfM;
      
df2_dT2=  - EM/(TM*h) + EM/TM^2*(T2-T1)/(2*h) - EM*HM/TM*(U2-U1)/(2*h*UM) ...
          + dE2_dT2/h + dE2_dT2*(T2-T1)/(2*TM*h) + dE2_dT2*(1+DM/TM)*(U2-U1)/(2*h*UM) ...
          + dE2_dT2/(2*TM)*CfM - dCD2_dT2/2 + EM/(2*TM)*dCf2_dT2 - EM/(2*TM^2)*CfM;

df2_dD1=  EM/TM*(U2-U1)/(2*h*UM) - dCD1_dD1/2 + EM/(2*TM)*dCf1_dD1 ...
         +dE1_dD1/TM*(U2-U1)/(2*h*UM) + dE1_dD1/(2*TM)*Cf1 ;
df2_dD2=  EM/TM*(U2-U1)/(2*h*UM) - dCD2_dD2/2 + EM/(2*TM)*dCf2_dD2 ...
         +dE2_dD2/TM*(U2-U1)/(2*h*UM) + dE2_dD2/(2*TM)*Cf2 ;

df2_dU1= - EM*(1+DM/TM)/(h*UM) - EM*(1+DM/TM)*(U2-U1)/(2*h*UM^2) ...
         - dCD1_dU1/2 +  EM/(2*TM)*dCf1_dU1;
df2_dU2=   EM*(1+DM/TM)/(h*UM) - EM*(1+DM/TM)*(U2-U1)/(2*h*UM^2) ...
         - dCD2_dU2/2 +  EM/(2*TM)*dCf2_dU2;

%  Jacobi-matrix point 1 and 2
%------------------------------------

% J1= [df1_dT1 , df1_dD1, df1_dU1 ;...
%      df2_dT1 , df2_dD1, df2_dU1 ];
%  
% J2= [df1_dT2 , df1_dD2, df1_dU2 ;...
%      df2_dT2 , df2_dD2, df2_dU2 ];

%  switch to m as primary variable
%------------------------------------


dU1_dm= D(1,:);
dU2_dm= D(2,:);
dD1_dm= 1/U1 - m1/U1^2 *dU1_dm; 
dD2_dm= 1/U2 - m2/U2^2 *dU2_dm; 

df1_dm=   df1_dD1*dD1_dm + df1_dD2*dD2_dm ...
        + df1_dU1*dU1_dm + df1_dU2*dU2_dm;
    
df2_dm=   df2_dD1*dD1_dm + df2_dD2*dD2_dm ...
        + df2_dU1*dU1_dm + df2_dU2*dU2_dm;
    
JT= [df1_dT1 , df1_dT2; df2_dT1, df2_dT2];

JM= [df1_dm ; df2_dm ];   
    
 
end

