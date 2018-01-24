function [ f1,f2,df_dT,df_dD,params, df_dU, df_dh] = SingleJacobiLam( DI,T,U,Vb,h,LTE,diffU,IsWake)
%SINGLEJACOBILAM    calculates the funtcionvalue and the derivates for Newton method
%                   of the momentum EQ and the shape parameter EQ in laminar case for known "1" values
%       Variables:  DI: displacementthicknes
%                   T : momentum thicknes
%                   U : tangential velocity at boundary edge
%                   Vb: normal velocity (uniform blowing)
%                   h : Vector with discretisation stepsize
%                   LTE: Trailing edge panel length
%                   diffU:  true returns also derivates of U
%                   IsWake: true if evaluation on wake


nu=evalin('base','nu');
N=length(T); 


if IsWake
    HWM = 0.5*(LTE(1:end-1)./T(1:end-1) + LTE(2:end)./T(2:end));  % Wake gap added to DI
    dHWM_dT2=-0.5*LTE(2:end)./T(2:end);
else
    HWM=zeros(size(h));
    dHWM_dT2=zeros(size(h));
end

% midpoint approximations
%----------------------------------------------
TM =0.5*( T(1:end-1) + T(2:end) ) ;
%DM =0.5*(DI(1:end-1) + DI(2:end)) ;
UM =0.5*( U(1:end-1) + U(2:end) ) ;
VbM=0.5*(Vb(1:end-1) + Vb(2:end)) ;



% help quantities H12 and ReT + derivates
%----------------------------------------------
H=DI./T; %shape parameter H1
H(H<1.05)=1.05;
HM=0.5*( H(1:end-1) + H(2:end) ); % midpoint value approximation
dH_dT   =-H ./T;
dH_dD   = 1./T;

Ret= U.*T/nu; 
RetM=0.5*( Ret(1:end-1) + Ret(2:end) ); % midpoint value approximation
dRet_dT = U/nu;
%dRet_dU = T/nu;



% shape parameter H32 + derivates
%----------------------------------------------
[ HS, dHS_dH ]=H32lam(H);
dHS_dT  = dHS_dH.*dH_dT;
dHS_dD  = dHS_dH.*dH_dD;


% friction coefficient cf/2
%----------------------------------------------
[ Cf, dCf_dH, dCf_dRet ]    = CF2lam( H, Ret);
[ CfM, dCfM_dH, dCfM_dRet ] = CF2lam( HM, RetM);
% special midpoint approximation
CfM=CfM/2 + (Cf(2:end) + Cf(1:end-1))/4;

dCf_dT= dCf_dH.*dH_dT + dCf_dRet.*dRet_dT;
dCf_dD= dCf_dH.*dH_dD;
%dCf_dU=               + dCf_dRet.*dRet_dU;


dCfM_dT2= 0.25*dCfM_dH.*dH_dT(2:end  ) + 0.25*dCfM_dRet.*dRet_dT(2:end  ) + 0.25*dCf_dT(2:end  );
%dCfM_dT1= 0.25*dCfM_dH.*dH_dT(1:end-1) + 0.25*dCfM_dRet.*dRet_dT(1:end-1) + 0.25*dCf_dT(1:end-1);

dCfM_dD2= 0.25*dCfM_dH.*dH_dD(2:end  ) + 0.25*dCf_dD(2:end  );
%dCfM_dD1= 0.25*dCfM_dH.*dH_dD(1:end-1)  + 0.25*dCf_dD(1:end-1);



% Dissipation coefficient 2 CD
%----------------------------------------------
[ CD, dCD_dH, dCD_dRet ]= CD2lam(H,Ret);

dCD_dT = dCD_dH.*dH_dT + dCD_dRet.*dRet_dT ;
dCD_dD = dCD_dH.*dH_dD ;
% dCD_dU =               + dCD_dRet.*dRet_dU;



% derivates of U if needed 
if diffU
    dRet_dU= T/nu; 
    dCf_dU = dCf_dRet.*dRet_dU;
    dCfM_dU2= 0.25*dCfM_dRet.*dRet_dU(2:end) + 0.25*dCf_dU(2:end) ;
    dCD_dU = dCD_dRet.*dRet_dU;
end



% use log(phi2/phi1) instead of (phi2-phi1)/phiM
%-------------------------------------------------------------------------
Ulog  = log( U(2:end)./U(1:end-1) ); % lin of U2-U1/UM
Tlog  = log( T(2:end)./T(1:end-1) );
HSlog = log(HS(2:end)./HS(1:end-1));

%--------------------------------------------------------------------------------
% Set up the Equations and the derivates for the Newton System
% momentum Equation divided by T and multiplied with h
% shape Parameter Equation divided by H32 and multiplied with h
%--------------------------------------------------------------------------------


% momentum equation -> right hand side of Newton System
f1= Tlog + (HM + 2 + HWM).*Ulog - h.*(CfM + VbM./UM)./TM;


% ------- derivates of momentum equation ------
% momentum thickness T
df1_dT2 = 1./T(2:end) + (dH_dT(2:end)/2 + dHWM_dT2).*Ulog - h.*dCfM_dT2./TM + 0.5*h.*(CfM + VbM./UM)./TM.^2;
% displacement thickness DI
df1_dD2 =             + (dH_dD(2:end)/2).*Ulog - h.*dCfM_dD2./TM ;
% Boundary layer edge velocity U
if diffU
    df1_dU2 =  (HM + 2 + HWM)./U(2:end) - h.*dCfM_dU2./TM + 0.5*h.*VbM./(UM.^2.*TM );
end
% panel length (only for Transition panel)
df1_dh=(CfM + VbM./UM)./TM;


%-------------------------------------------------------------------------


% Upwinding for shape parameter equation -> the more the value of H12 at 1 differs from the value at 2, the more upwinding
% upw= 0.5 -> central value approximation, upw = 1 -> backward approximation (2 values)
[upw, ~ ,dupw_dH2]=Upwinding(H(1:end-1),H(2:end),false);
dupw_dT2=dupw_dH2.*dH_dT(2:end  );
dupw_dD2=dupw_dH2.*dH_dD(2:end  );

I=ones(N-1,1);
HSM=(I-upw).*HS(1:end-1) + upw.*HS(2:end);
CfM=(I-upw).*Cf(1:end-1) + upw.*Cf(2:end);
CDM=(I-upw).*CD(1:end-1) + upw.*CD(2:end);

dHSM_dT2= upw    .*dHS_dT(2:end  ) + (HS(2:end)-HS(1:end-1)).*dupw_dT2;
dHSM_dD2= upw    .*dHS_dD(2:end  ) + (HS(2:end)-HS(1:end-1)).*dupw_dD2;
dCfM_dT2= upw    .*dCf_dT(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT2;
dCfM_dD2= upw    .*dCf_dD(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD2;
dCDM_dT2= upw    .*dCD_dT(2:end  ) + (CD(2:end)-CD(1:end-1)).*dupw_dT2;
dCDM_dD2= upw    .*dCD_dD(2:end  ) + (CD(2:end)-CD(1:end-1)).*dupw_dD2;

if diffU
    dCfM_dU2= upw    .*dCf_dU(2:end  );
    dCDM_dU2= upw    .*dCD_dU(2:end  );
end


% shape parameter equation -> right hand side of Newton System
f2= HSlog + (1-HM-HWM).*Ulog + h.*(CfM - CDM + VbM./UM.*(1-1./HSM) )./TM;


% derivates of shape parameter equation
%-------------------------------------------------------------------------

% momentum thickness T
df2_dT2 = 1./HS(2:end).*dHS_dT(2:end) - (dH_dT(2:end)/2 + dHWM_dT2).*Ulog ...
     + h.*(dCfM_dT2 - dCDM_dT2 )./TM - 0.5*h.*(CfM - CDM + VbM./UM.*(1-1./HSM) )./TM.^2 + VbM./(UM.*TM.*HSM.^2).*dHSM_dT2; 
% displacement thickness DI
df2_dD2 = 1./HS(2:end).*dHS_dD(2:end) - (dH_dD(2:end)/2).*Ulog + h.*(dCfM_dD2 - dCDM_dD2)./TM + VbM./(UM.*TM.*HSM.^2).*dHSM_dD2 ; 
% Boundary layer edge velocity U
if diffU
    df2_dU2 =   (1-HM-HWM)./U(2:end) + h.*(dCfM_dU2 - dCDM_dU2)./TM - 0.5*h.*VbM.*(1-1./HSM)./(UM.^2.*TM )  ; 
    df_dU=[df1_dU2; df2_dU2]; 
end
% panel length (only for Transition panel)
df2_dh=(CfM - CDM + VbM./UM)./TM;


df_dD=[df1_dD2; df2_dD2]; 
df_dT=[df1_dT2; df2_dT2]; 
df_dh=[df1_dh; df2_dh]; 

params= [CfM, Cf(1), CDM, CD(1), HS(1), Ret(1) ; 0, Cf(2), 0, CD(2),HS(2), Ret(2) ];
end

