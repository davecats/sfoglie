function [f1,f2, J ] = InitialNodeSys(T,D,U,Vb,h,nu)
%INITIALNODESYS sets up a Newton System for first Boundary node. The
%equations are those of the intervall between Stagnation point and
%first boundary layer node setting all quantities to 0 at the stagnation
%point and using a backward difference scheme.

 
H=D/T; 
dH_dT   =-H /T;
dH_dD   = 1/T;

Ret= T*U/nu;
dRet_dT = U/nu;
dRet_dU = T/nu;

[ HS, dHS_dH ]=H32lam(H);
dHS_dT  = dHS_dH*dH_dT;
dHS_dD  = dHS_dH*dH_dD;


[ Cf, dCf_dH, dCf_dRet ]    = CF2lam( H, Ret);
dCf_dT= dCf_dH*dH_dT + dCf_dRet*dRet_dT;
dCf_dD= dCf_dH*dH_dD;
dCf_dU=               + dCf_dRet*dRet_dU;


[ CD, dCD_dH, dCD_dRet ]= CD2lam(H,Ret);
dCD_dT = dCD_dH*dH_dT + dCD_dRet*dRet_dT ;
dCD_dD = dCD_dH*dH_dD ;
dCD_dU =               + dCD_dRet*dRet_dU;

f1=  (H + 2 ) - h.*(Cf + Vb./U)./T;

% momentum thickness T
df1_dT =    dH_dT - h*dCf_dT/T  + h*(Cf + Vb/U)/T^2;

% displacement thickness DI
df1_dD =    dH_dD - h*dCf_dD/T ;

% Boundary layer edge velocity U
df1_dU =   - h*dCf_dU/T + h*Vb/(U^2*T );


df1_ds =   - (Cf + Vb./U)./T;


f2= (1-H) + h*(Cf - CD + Vb/U*(1-1/HS) )/T;

% derivates of shape parameter equation
%-------------------------------------------------------------------------

% momentum thickness T
df2_dT =  - dH_dT + h*(dCf_dT - dCD_dT )/T - h*(Cf - CD + Vb/U*(1-1/HS) )/T^2 + h*Vb/(U*T*HS^2).*dHS_dT; 

% displacement thickness DI
df2_dD =  - dH_dD + h*(dCf_dD - dCD_dD)/T + h*Vb/(U*T*HS^2)*dHS_dD ; 

% Boundary layer edge velocity U
df2_dU =    h*(dCf_dU - dCD_dU)/T - h*Vb*(1-1/HS)/(U^2*T)  ; 


df2_ds =  (Cf - CD + Vb/U*(1-1/HS) )/T;


J=[ df1_dT , df1_dD, df1_dU , df1_ds; df2_dT , df2_dD, df2_dU , df2_ds];


end

