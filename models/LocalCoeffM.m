function [ JT, JM,f1, f2 ] = LocalCoeffM(i,h,T,Uinv,m,D)
%LOCALNEWTONCOEFF   calculates the contribution of each panel to the global Newton System
%                   panel with delta_2 and mass defect m as primary variables
%                   i: panel index
%                   h: panel length
%                   T: Vector with momentum thickness of start and end point
%                   Uinv: Vector with inviscid velocity of start and end point
%                   m: Vector of all mass defects
%                   D: Coefficients of mass Defekt -> U = Uinv + Dm 
%                           -> i-th and i+1-th line

nu=evalin('base','nu');
%nu=1.53e-5;

% momentum thickness
T1=T(1);
T2=T(1);
% bouandary layer edge velocity
U1=D(1,:)*m +Uinv(1);
U2=D(2,:)*m +Uinv(2);

% displacement thickness
D1=m(i)/U1;
D2=m(i+1)/U2;


% panel midpoint values approximation
TM=(T1+T2)/2; % momentum thickness
UM=(U1+U2)/2; % Velocity
DM=(D1+D2)/2; % displacement thickness
mM=(m(i+1)+m(i))/2; % mass defect
H12=DM/TM;    % shape parameter  
Ret=UM*TM/nu; % mom. thickn. Reynolds-number

% derivates of D,U,H12 in respect to T and m
%dU_dm= (D(i,i)+D(i+1,i+1) )/2;
dU_dm= (D(1,:)+D(2,:) )/2;
dH_dT=-H12/TM;
dD_dm=1/UM - mM/UM^2 *dU_dm; 
dH_dm= dD_dm/TM;

dRet_dT= UM/nu;
dRet_dm= TM/nu * dU_dm;



%  model for Cf/2
%------------------------------------
[Cf2, dCf2_dH, dCf2_dRet]=CF2lam(H12, Ret);

% derivates of primary variables in respect to T and m
dCf2_dT = dCf2_dH*dH_dT + dCf2_dRet*dRet_dT  ;
dCf2_dm = dCf2_dH*dH_dm + dCf2_dRet*dRet_dm  ;


%  model for cinetic energy thickness (ET)
%------------------------------------

% mean values

[EM,dEM_dH]=delta3lam(TM,H12);

% derivates of primary variables
dEM_dT= dEM_dH*dH_dT + EM ;
dEM_dm= dEM_dH*dH_dm;

% difference values
[E1,dE1_dH]=delta3lam(T1,D1/T1);
dE1_dT= dE1_dH*dH_dT;
dE1_dm= dE1_dH*dH_dm;

[E2,dE2_dH]=delta3lam(T2,D2/T2);
dE2_dT= dE2_dH*dH_dT;
dE2_dm= dE2_dH*dH_dm;

dH32_dT= dEM_dT/TM - EM/TM^2;
dH32_dm= 1/TM*dEM_dm;

%  model for dissipation integral 2*CD
%------------------------------------

[CD2, dCD2_dH,dCD2_dRet]=CD2lam(H12,Ret,EM/TM);

% derivates of primary variables
dCD2_dT=dCD2_dH*dH_dT + dCD2_dRet*dRet_dT + CD2*dH32_dT;
dCD2_dm=dCD2_dH*dH_dm + dCD2_dRet*dRet_dm + CD2*dH32_dm;






%  momentum Equation
%------------------------------------
f1=(T2-T1)/h+(2*TM+DM)*1/UM*(U2-U1)/h-Cf2;

% derivates of momentum Eq in respect to primary variables
df1_dT1=-1/h + (U2-U1)/(h*UM) - dCf2_dT /2;
df1_dT2= 1/h + (U2-U1)/(h*UM) - dCf2_dT /2;

dU1_dm= D(1,:);
dU2_dm= D(2,:);

df1_dm = dD_dm*(U2-U1)/(h*UM) - dU_dm*(2*TM+DM)*(U2-U1)/(h*UM^2) ...
            +(dU2_dm-dU1_dm)*(2*TM+DM)/UM - dCf2_dm;

%df1_dm1 = dD_dm*(U2-U1)/(2*h*UM) - dU_dm*(2*TM+DM)*(U2-U1)/(h*UM^2) ...
%            +(dU2_dm-dU1_dm)*(2*TM+DM)/UM - dCf2_dm;

%  shape parameter Equation
%------------------------------------

f2= (E2-E1)/h - EM/TM* (T2-T1)/h + EM*(1+H12)*1/UM*(U2-U1)/h -CD2 + EM/TM*Cf2;

% derivates of shape parameter Eq in respect to primary variables
df2_dT1= (dE2_dT-dE1_dT)/(h) -dEM_dT/TM*(T2-T1)/h - dEM_dT*(1+H12)*1/UM*(U2-U1)/h ...
          + dEM_dT/TM*Cf2 + EM/TM^2*(T2-T1)/h + EM/(h*TM) - EM/TM^2*Cf2 ...
          + dH_dT*EM/UM*(U2-U1)/h - dCD2_dT + EM/TM*dCf2_dT;

df2_dT2= (dE2_dT-dE1_dT)/(h) -dEM_dT/TM*(T2-T1)/h - dEM_dT*(1+H12)*1/UM*(U2-U1)/h ...
          + dEM_dT/TM*Cf2 + EM/TM^2*(T2-T1)/h - EM/(h*TM) - EM/TM^2*Cf2 ...
          + dH_dT*EM/UM*(U2-U1)/h - dCD2_dT + EM/TM*dCf2_dT;

      

df2_dm= (dE2_dm-dE1_dm)/h - dEM_dm/TM* (T2-T1)/h + dEM_dm*(1+DM/TM)*1/UM*(U2-U1)/h ...
         + dEM_dm/TM*Cf2 + dH_dm*EM*1/UM*(U2-U1)/h - dU_dm*EM*(1+H12)*1/UM^2*(U2-U1)/h ...
         + EM*(1+H12)*1/UM*(dU2_dm-dU1_dm)/h -dCD2_dm + EM/TM*dCf2_dm ;
 
 %{    
 % df2_dT1=    1/(2*h)*ET*(T2-T1)/TM^2 + ET/(TM*h) - 1/(2*h)*ET*DM*(U2-U1)/(TM^2*UM) ...
%           + (dET2_dT-dET1_dT)/h + dET_dT/2*(T2-T1)/(TM*h) + dET_dT/2*(1+DM/TM)*(U2-U1)/(UM*h)...
%           +  dET_dT/(2*TM)*Cf2 - dCD2_dT - ET/TM^2*Cf2 + ET/TM*dCf_dT;
%       
%       
% df2_dT2=  +  1/(2*h)*ET*(T2-T1)/TM^2 - ET/(TM*h) - 1/(2*h)*ET*DM*(U2-U1)/(TM^2*UM) ...
%           + (dET2_dT-dET1_dT)/h + dET_dT/2*(T2-T1)/(TM*h) + dET_dT/2*(1+DM/TM)*(U2-U1)/(UM*h)...
%           +  dET_dT/(2*TM)*Cf2 - dCD2_dT - ET/TM^2*Cf2 + ET/TM*dCf_dT;   
%}     
     
 %  Jacobi-matrix point 1 and 2
%------------------------------------

JT= [df1_dT1 , df1_dT2; df2_dT1, df2_dT2];

%JM= [df1_dm1, df1_dm2 ; df2_dm1, df2_dm2 ];
JM= [df1_dm ; df2_dm ];




end

