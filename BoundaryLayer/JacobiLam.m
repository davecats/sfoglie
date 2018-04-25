function [ f1,f2,der] = JacobiLam( T,U,Vb,H,Ret,h,s1,nu,pressureTerm)
%JACOBILAM          calculates the funtcion value and the derivates for Newton method
%                   of the momentum EQ and the shape parameter EQ in
%                   laminar case for start and end point of each intervall
%       Variables:  T : momentum thicknes
%                   U : tangential velocity at boundary edge
%                   Vb: normal velocity (uniform blowing)
%                   H : shape parameter H12
%                   h : Vector with discretisation stepsize
%                   pressureTerm: correction of pressure only in momentum eq

N=length(T); 
if nargin==8
    pressureTerm=zeros(size(Vb));
end
pTM=0.5*(pressureTerm(1:end-1) + pressureTerm(2:end)) ;

s2=s1+h;
sM=0.5*(s2+s1);
T1=T(1:end-1);
T2=T(2:end);

% midpoint approximations
%----------------------------------------------
TM =0.5*( T(1:end-1) + T(2:end) ) ;
UM =0.5*( U(1:end-1) + U(2:end) ) ;
VbM=0.5*(Vb(1:end-1) + Vb(2:end)) ;



% help quantities H12 and ReT + derivates
%----------------------------------------------
HM=0.5*( H(1:end-1) + H(2:end) ); % midpoint value approximation
dH_dT   =-H ./T;
dH_dD   = 1 ./T;

RetM=0.5*( Ret(1:end-1) + Ret(2:end) ); % midpoint value approximation
dRet_dT = U/nu;
dRet_dU = T/nu;



% shape parameter H32 + derivates
%----------------------------------------------
[ HS, dHS_dH ]=H32lam(H);
dHS_dT  = dHS_dH.*dH_dT;
dHS_dD  = dHS_dH.*dH_dD;

% friction coefficient cf/2
%----------------------------------------------
[ Cf, dCf_dH, dCf_dRet ]    = CF2lam( H, Ret);
[ CfM, dCfM_dH, dCfM_dRet ] = CF2lam( HM, RetM);


dCf_dT= dCf_dH.*dH_dT + dCf_dRet.*dRet_dT;
dCf_dD= dCf_dH.*dH_dD;
dCf_dU=               + dCf_dRet.*dRet_dU;


dCfM_dT2= 0.5*dCfM_dH.*dH_dT(2:end  ) + 0.5*dCfM_dRet.*dRet_dT(2:end  );
dCfM_dT1= 0.5*dCfM_dH.*dH_dT(1:end-1) + 0.5*dCfM_dRet.*dRet_dT(1:end-1);

dCfM_dD2= 0.5*dCfM_dH.*dH_dD(2:end  )  ;
dCfM_dD1= 0.5*dCfM_dH.*dH_dD(1:end-1)  ;

dCfM_dU2= 0.5*dCfM_dRet.*dRet_dU(2:end)   ;
dCfM_dU1= 0.5*dCfM_dRet.*dRet_dU(1:end-1) ;

% incorporate s/T already in CF value -> better for initial steps to prevent oscillations
sT1=s1./T1;
sT2=s2./T2;
sTM=sM./TM;
CFX=0.5*CfM.*sTM + 0.25*(Cf(2:end).*sT2 + Cf(1:end-1).*sT1);
dCFX_dT2= 0.5*dCfM_dT2.*sTM + 0.25*(dCf_dT(2:end  ).*sT2) - 0.25*CfM.*sTM./TM - 0.25*Cf(2:end)  .*sT2./T2;
dCFX_dT1= 0.5*dCfM_dT1.*sTM + 0.25*(dCf_dT(1:end-1).*sT1) - 0.25*CfM.*sTM./TM - 0.25*Cf(1:end-1).*sT1./T1;
dCFX_dD2= 0.5*dCfM_dD2.*sTM + 0.25*dCf_dD(2:end  ).*sT2;
dCFX_dD1= 0.5*dCfM_dD1.*sTM + 0.25*dCf_dD(1:end-1).*sT1;
dCFX_dU2= 0.5*dCfM_dU2.*sTM + 0.25*dCf_dU(2:end  ).*sT2;
dCFX_dU1= 0.5*dCfM_dU1.*sTM + 0.25*dCf_dU(1:end-1).*sT1;

dCFX_ds2= 0.25*CfM./TM + 0.25*Cf(2:end  )./T2;
dCFX_ds1= 0.25*CfM./TM + 0.25*Cf(1:end-1)./T1;


% Dissipation coefficient 2 CD
%----------------------------------------------
[ CD, dCD_dH, dCD_dRet ]= CD2lam(H,Ret);

dCD_dT = dCD_dH.*dH_dT + dCD_dRet.*dRet_dT ;
dCD_dD = dCD_dH.*dH_dD ;
dCD_dU =               + dCD_dRet.*dRet_dU;


% use log(phi2/phi1) instead of (phi2-phi1)/phiM
%-------------------------------------------------------------------------
Ulog  = log( U(2:end)./U(1:end-1) ); % lin of U2-U1/UM
Tlog  = log( T(2:end)./T(1:end-1) );
HSlog = log(HS(2:end)./HS(1:end-1));
slog  = log(s2./s1); % substitute h= log(s2/s1) * sM


%--------------------------------------------------------------------------------
% Set up the Equations and the derivates for the Newton System
% momentum Equation divided by T and multiplied with h
% shape Parameter Equation divided by H32 and multiplied with h
%--------------------------------------------------------------------------------


% momentum equation -> right hand side of Newton System
f1= Tlog + (HM + 2 ).*Ulog - slog.*CFX - h.*VbM./(UM.*TM) - h.*pTM./(UM.^2.*TM);


% ------- derivates of momentum equation ------
% momentum thickness T
df1_dT2 =  1./T2 + 0.5*dH_dT(2:end  ).*Ulog - slog.*dCFX_dT2 + 0.5*h.* VbM./(UM.*TM.^2) + 0.5*h.* pTM./(UM.^2.*TM.^2);
df1_dT1 = -1./T1 + 0.5*dH_dT(1:end-1).*Ulog - slog.*dCFX_dT1 + 0.5*h.* VbM./(UM.*TM.^2) + 0.5*h.* pTM./(UM.^2.*TM.^2);

% displacement thickness DI
df1_dD2 =                  0.5*dH_dD(2:end  ).*Ulog - slog.*dCFX_dD2 ;
df1_dD1 =                  0.5*dH_dD(1:end-1).*Ulog - slog.*dCFX_dD1 ;

% Boundary layer edge velocity U
df1_dU2 =   (HM + 2 )./U(2:end  ) - slog.*dCFX_dU2 + 0.5*h.*VbM./(UM.^2.*TM ) + h.*pTM./(UM.^3.*TM );
df1_dU1 =  -(HM + 2 )./U(1:end-1) - slog.*dCFX_dU1 + 0.5*h.*VbM./(UM.^2.*TM ) + h.*pTM./(UM.^3.*TM );
% arc length s
df1_ds2=-CFX./s2 - slog.*dCFX_ds2 - VbM./(UM.*TM) - pTM./(UM.^2.*TM);
df1_ds1= CFX./s1 - slog.*dCFX_ds1 + VbM./(UM.*TM) + pTM./(UM.^2.*TM);


%-------------------------------------------------------------------------


% Upwinding for shape parameter equation -> the more the value of H12 at 1 differs from the value at 2, the more upwinding
% upw= 0.5 -> central value approximation, upw = 1 -> backward approximation (2 values)
[upw, dupw_dH1 ,dupw_dH2]=Upwinding(H(1:end-1),H(2:end),false);
dupw_dT2=dupw_dH2.*dH_dT(2:end  );
dupw_dD2=dupw_dH2.*dH_dD(2:end  );
dupw_dT1=dupw_dH1.*dH_dT(1:end-1);
dupw_dD1=dupw_dH1.*dH_dD(1:end-1);


I=ones(N-1,1);
HSM=(I-upw).*HS(1:end-1) + upw.*HS(2:end);
% CfM=(I-upw).*Cf(1:end-1) + upw.*Cf(2:end);
% CDM=(I-upw).*CD(1:end-1) + upw.*CD(2:end);
CFX=(I-upw).*Cf(1:end-1).*sT1 + upw.*Cf(2:end).*sT2;
CDX=(I-upw).*CD(1:end-1).*sT1 + upw.*CD(2:end).*sT2;


dHSM_dT2= upw    .*dHS_dT(2:end  ) + (HS(2:end)-HS(1:end-1)).*dupw_dT2;
dHSM_dD2= upw    .*dHS_dD(2:end  ) + (HS(2:end)-HS(1:end-1)).*dupw_dD2;
% dCfM_dT2= upw    .*dCf_dT(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT2;
% dCfM_dD2= upw    .*dCf_dD(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD2;
% dCDM_dT2= upw    .*dCD_dT(2:end  ) + (CD(2:end)-CD(1:end-1)).*dupw_dT2;
% dCDM_dD2= upw    .*dCD_dD(2:end  ) + (CD(2:end)-CD(1:end-1)).*dupw_dD2;
% dCfM_dU2= upw    .*dCf_dU(2:end  );
% dCDM_dU2= upw    .*dCD_dU(2:end  );

dCFX_dT2= upw    .*(dCf_dT(2:end  ).*sT2 - Cf(2:end  ).*sT2./T2) + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dT2;
dCFX_dD2= upw    .*(dCf_dD(2:end  ).*sT2)                        + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dD2;    
dCDX_dT2= upw    .*(dCD_dT(2:end  ).*sT2 - CD(2:end  ).*sT2./T2) + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dT2;
dCDX_dD2= upw    .*(dCD_dD(2:end  ).*sT2)                        + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dD2; 
dCFX_dU2= upw    .*(dCf_dU(2:end  ).*sT2);
dCDX_dU2= upw    .*(dCD_dU(2:end  ).*sT2);
dCFX_ds2= upw    .*(Cf(2:end  )./T2);
dCDX_ds2= upw    .*(CD(2:end  )./T2);

dHSM_dT1= (I-upw).*dHS_dT(1:end-1) + (HS(2:end)-HS(1:end-1)).*dupw_dT1;
dHSM_dD1= (I-upw).*dHS_dD(1:end-1) + (HS(2:end)-HS(1:end-1)).*dupw_dD1;
% dCfM_dT1= (I-upw).*dCf_dT(1:end-1) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT1;
% dCfM_dD1= (I-upw).*dCf_dD(1:end-1) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD1;
% dCDM_dT1= (I-upw).*dCD_dT(1:end-1) + (CD(2:end)-CD(1:end-1)).*dupw_dT1;
% dCDM_dD1= (I-upw).*dCD_dD(1:end-1) + (CD(2:end)-CD(1:end-1)).*dupw_dD1;
% dCfM_dU1= (I-upw).*dCf_dU(1:end-1);
% dCDM_dU1= (I-upw).*dCD_dU(1:end-1);
dCFX_dT1= (I-upw).*(dCf_dT(1:end-1).*sT1 - Cf(1:end-1).*sT1./T1) + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dT1;
dCFX_dD1= (I-upw).*(dCf_dD(1:end-1).*sT1)                        + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dD1;    
dCDX_dT1= (I-upw).*(dCD_dT(1:end-1).*sT1 - CD(1:end-1).*sT1./T1) + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dT1;
dCDX_dD1= (I-upw).*(dCD_dD(1:end-1).*sT1)                        + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dD1; 
dCFX_dU1= (I-upw).*(dCf_dU(1:end-1).*sT1);
dCDX_dU1= (I-upw).*(dCD_dU(1:end-1).*sT1);
dCFX_ds1= (I-upw).*(Cf(1:end-1)./T1);
dCDX_ds1= (I-upw).*(CD(1:end-1)./T1);



% shape parameter equation -> right hand side of Newton System
f2= HSlog + (1-HM).*Ulog + slog.*(CFX - CDX) + h.*VbM./UM.*(1-1./HSM) ./TM;


% derivates of shape parameter equation
%-------------------------------------------------------------------------

% momentum thickness T
df2_dT2 =  1./HS(2:end  ).*dHS_dT(2:end  ) - 0.5*dH_dT(2:end  ).*Ulog ...
            + slog.*(dCFX_dT2 - dCDX_dT2 ) - 0.5*h.*VbM./UM.*(1-1./HSM)./TM.^2 + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dT2; 
 
df2_dT1 = -1./HS(1:end-1).*dHS_dT(1:end-1) - 0.5*dH_dT(1:end-1).*Ulog ...
            + slog.*(dCFX_dT1 - dCDX_dT1 ) - 0.5*h.*VbM./UM.*(1-1./HSM)./TM.^2 + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dT1; 
 
 
% displacement thickness DI
df2_dD2 = 1./HS(2:end  ).*dHS_dD(2:end  ) - 0.5*dH_dD(2:end  ).*Ulog + slog.*(dCFX_dD2 - dCDX_dD2) + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dD2 ; 
df2_dD1 =-1./HS(1:end-1).*dHS_dD(1:end-1) - 0.5*dH_dD(1:end-1).*Ulog + slog.*(dCFX_dD1 - dCDX_dD1) + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dD1 ; 

% Boundary layer edge velocity U
df2_dU2 =   (1-HM)./U(2:end  ) + slog.*(dCFX_dU2 - dCDX_dU2) - 0.5*h.*VbM.*(1-1./HSM)./(UM.^2.*TM );
df2_dU1 =  -(1-HM)./U(1:end-1) + slog.*(dCFX_dU1 - dCDX_dU1) - 0.5*h.*VbM.*(1-1./HSM)./(UM.^2.*TM );

% arc length s
df2_ds2= 1./s2.*(CFX - CDX) + slog.*(dCFX_ds2 - dCDX_ds2) + VbM./UM.*(1-1./HSM) ./TM;
df2_ds1=-1./s1.*(CFX - CDX) + slog.*(dCFX_ds1 - dCDX_ds1) - VbM./UM.*(1-1./HSM) ./TM;

der.df1_dD=[df1_dD1, df1_dD2]; 
der.df1_dT=[df1_dT1, df1_dT2]; 
der.df1_dU=[df1_dU1, df1_dU2];  
der.df2_dD=[df2_dD1, df2_dD2]; 
der.df2_dT=[df2_dT1, df2_dT2]; 
der.df2_dU=[df2_dU1, df2_dU2];  

der.df1_ds=[df1_ds1, df1_ds2];
der.df2_ds=[df2_ds1, df2_ds2];
end

