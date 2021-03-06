function [ f1,f2,f3,der ] = JacobiTurb( DI, T,Ct,U,Vb,H,Ret,h,s1,IsWake,Wgap,nu,pressureTerm)
%SINGLEJACOBITURB   calculates the funtcion value and the derivates of
%                   turbulent part  for the Newton method.
%                   f1: momentum EQ 
%                   f2: shape parameter EQ 
%                   f3: lag entrainement EQ (Ctau)    
%       Variables:  DI: displacement thickness
%                   T : momentum thickness
%                   Ct: max shear stress coefficient Ctau
%                   U : tangential velocity at boundary edge
%                   Vb: normal velocity (uniform blowing)
%                   h : Vector with discretisation stepsize
%                   IsWake: true if evaluation on wake



N=length(T); 
dimPl=size(T);
dimMi=size(h);

if nargin==12
    pressureTerm=zeros(size(Vb));
end
pTM=0.5*(pressureTerm(1:end-1) + pressureTerm(2:end)) ;
    
if ~IsWake ; Wgap=zeros(size(U)); end

s2=s1+h;
sM=(s2+s1)./2;
T1=T(1:end-1);
T2=T(2:end);

% substract TE Gap for wake values
DI=DI-Wgap;
HW=Wgap./T;
HWM=0.5*( HW(1:end-1) + HW(2:end) );
dHWM_dT=-0.5*HW./T;


% midpoint approximations
%----------------------------------------------
TM =0.5*( T(1:end-1) + T(2:end) ) ;
DM =0.5*(DI(1:end-1) + DI(2:end)) ;
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
[ HS, dHS_dH, dHS_dRet  ]=H32turb(H,Ret);

dHS_dT  = dHS_dH.*dH_dT + dHS_dRet.*dRet_dT;
dHS_dD  = dHS_dH.*dH_dD;
dHS_dU  =                dHS_dRet.*dRet_dU;

% friction coefficient cf/2
%----------------------------------------------
if IsWake
    % no friction coeff at wake
    Cf=zeros(dimPl); dCf_dH=zeros(dimPl); dCf_dRet=zeros(dimPl);
    CfM=zeros(dimMi);dCfM_dH=zeros(dimMi);dCfM_dRet =zeros(dimMi);
else
    [ Cf, dCf_dH, dCf_dRet ]    = CF2turb( H, Ret);
    [ CfM, dCfM_dH, dCfM_dRet ] = CF2turb( HM, RetM);
end


dCf_dT= dCf_dH.*dH_dT + dCf_dRet.*dRet_dT;
dCf_dD= dCf_dH.*dH_dD;
dCf_dU=                dCf_dRet.*dRet_dU;

dCfM_dT2= 0.5*dCfM_dH.*dH_dT(2:end  ) + 0.5*dCfM_dRet.*dRet_dT(2:end  ) ;
dCfM_dT1= 0.5*dCfM_dH.*dH_dT(1:end-1) + 0.5*dCfM_dRet.*dRet_dT(1:end-1) ;

dCfM_dD2= 0.5*dCfM_dH.*dH_dD(2:end  ) ;
dCfM_dD1= 0.5*dCfM_dH.*dH_dD(1:end-1) ;

dCfM_dU2= 0.5*dCfM_dRet.*dRet_dU(2:end  )  ;
dCfM_dU1= 0.5*dCfM_dRet.*dRet_dU(1:end-1)  ;

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

% Slip velocity Us
%----------------------------------------------
Us=0.5*HS.*( -1/3 + 1./(0.75*H) );
dUs_dHS=  Us./HS;
dUs_dH = -HS./(1.5*H.^2);

% filter high values
if IsWake
    ind=find(Us>0.99995);
    Us(ind)=0.99995; dUs_dHS(ind)=0; dUs_dH(ind)=0;
else
    ind=find(Us>0.95);
    Us(ind)=0.98; dUs_dHS(ind)=0; dUs_dH(ind)=0; 
end



USM=0.5*( Us(2:end)+Us(1:end-1) );

dUs_dT = dUs_dH.*dH_dT + dUs_dHS.*dHS_dT;
dUs_dD = dUs_dH.*dH_dD + dUs_dHS.*dHS_dD;
dUs_dU =                 dUs_dHS.*dHS_dU;

% Equilibrium Shear stress coefficient CE = sqrt(CtauEQ)
%----------------------------------------------

[ CE,dCE_dH,dCE_dRet,dCE_dHS,dCE_dUs ]=CtEQ( H,Ret,HS,Us,IsWake);
dCE_dT = dCE_dH.*dH_dT + dCE_dRet.*dRet_dT + dCE_dHS.*dHS_dT + dCE_dUs.*dUs_dT ;
dCE_dD = dCE_dH.*dH_dD +                   + dCE_dHS.*dHS_dD + dCE_dUs.*dUs_dD ;
dCE_dU =                 dCE_dRet.*dRet_dU + dCE_dHS.*dHS_dU + dCE_dUs.*dUs_dU ;


% Dissipation coefficient 2 CD
%----------------------------------------------
[ CD, dCD_dH,dCD_dRet, dCD_dHS,dCD_dUs, dCD_dCt] = CD2turb( H,Ret,HS,Us,Ct,IsWake);

dCD_dT = dCD_dH.*dH_dT + dCD_dRet.*dRet_dT + dCD_dHS.*dHS_dT + dCD_dUs.*dUs_dT  ;
dCD_dD = dCD_dH.*dH_dD +                   + dCD_dHS.*dHS_dD + dCD_dUs.*dUs_dD  ;
dCD_dU =                 dCD_dRet.*dRet_dU + dCD_dHS.*dHS_dU + dCD_dUs.*dUs_dU  ;



if IsWake % in case of wake 2 sides of boundary Layer -> double dissipation
   CD=2*CD;  
   dCD_dT=2*dCD_dT; 
   dCD_dD=2*dCD_dD; 
   dCD_dU=2*dCD_dU; 
   dCD_dCt=2*dCD_dCt;
end
    
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
f1= Tlog + (HM + HWM + 2 ).*Ulog - slog.*CFX - h.*VbM./(UM.*TM) - h.*pTM./(UM.^2.*TM);


% ------- derivates of momentum equation ------
% momentum thickness T
df1_dT2 =  1./T2 + (0.5*dH_dT(2:end  ) + dHWM_dT(2:end  )).*Ulog - slog.*dCFX_dT2 + 0.5*h.* VbM./(UM.*TM.^2) + 0.5*h.* pTM./(UM.^2.*TM.^2);
df1_dT1 = -1./T1 + (0.5*dH_dT(1:end-1) + dHWM_dT(1:end-1)).*Ulog - slog.*dCFX_dT1 + 0.5*h.* VbM./(UM.*TM.^2) + 0.5*h.* pTM./(UM.^2.*TM.^2);

% displacement thickness DI
df1_dD2 =           0.5*dH_dD(2:end  ).*Ulog                     - slog.*dCFX_dD2 ;
df1_dD1 =           0.5*dH_dD(1:end-1).*Ulog                     - slog.*dCFX_dD1 ;

% Boundary layer edge velocity U
df1_dU2 =   (HM + 2 + HWM)./U(2:end  ) - slog.*dCFX_dU2 + 0.5*h.*VbM./(UM.^2.*TM ) + h.*pTM./(UM.^3.*TM );
df1_dU1 =  -(HM + 2 + HWM)./U(1:end-1) - slog.*dCFX_dU1 + 0.5*h.*VbM./(UM.^2.*TM ) + h.*pTM./(UM.^3.*TM );

% arc length s
if IsWake
   df1_ds1=zeros(dimMi); df1_ds2=zeros(dimMi); 
else
   df1_ds2=-1./s2.*CFX - slog.*dCFX_ds2 - VbM./(UM.*TM) - pTM./(UM.^2.*TM);
   df1_ds1= 1./s1.*CFX - slog.*dCFX_ds1 + VbM./(UM.*TM) + pTM./(UM.^2.*TM);
end


%-------------------------------------------------------------------------


% Upwinding for shape parameter equation -> the more the value of H12 at 1 differs from the value at 2, the more upwinding
% upw= 0.5 -> central value approximation, upw = 1 -> backward approximation (2 values)
% Upwinding for shape parameter equation -> the more the value of H12 at 1 differs from the value at 2, the more upwinding
% upw= 0.5 -> central value approximation, upw = 1 -> backward approximation (2 values)
[upw, dupw_dH1 ,dupw_dH2]=Upwinding(H(1:end-1),H(2:end),IsWake);
dupw_dT2=dupw_dH2.*dH_dT(2:end  );
dupw_dD2=dupw_dH2.*dH_dD(2:end  );
dupw_dT1=dupw_dH1.*dH_dT(1:end-1);
dupw_dD1=dupw_dH1.*dH_dD(1:end-1);


I=ones(N-1,1);

HSM=(I-upw).*HS(1:end-1) + upw.*HS(2:end);

CFX=(I-upw).*Cf(1:end-1).*sT1 + upw.*Cf(2:end).*sT2;
CDX=(I-upw).*CD(1:end-1).*sT1 + upw.*CD(2:end).*sT2;

dHSM_dT2= upw    .*dHS_dT(2:end  ) + (HS(2:end)-HS(1:end-1)).*dupw_dT2;
dHSM_dD2= upw    .*dHS_dD(2:end  ) + (HS(2:end)-HS(1:end-1)).*dupw_dD2;
dHSM_dU2= upw    .*dHS_dU(2:end  );
dCFX_ds2= upw    .*(Cf(2:end  )./T2);
dCDX_ds2= upw    .*(CD(2:end  )./T2);

dHSM_dT1= (I-upw).*dHS_dT(1:end-1) + (HS(2:end)-HS(1:end-1)).*dupw_dT1;
dHSM_dD1= (I-upw).*dHS_dD(1:end-1) + (HS(2:end)-HS(1:end-1)).*dupw_dD1;
dHSM_dU1= (I-upw).*dHS_dU(1:end-1);
dCFX_ds1= (I-upw).*(Cf(1:end-1)./T1);
dCDX_ds1= (I-upw).*(CD(1:end-1)./T1);



dCFX_dT2= upw    .*(dCf_dT(2:end  ).*sT2 - Cf(2:end  ).*sT2./T2) + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dT2;
dCFX_dD2= upw    .*(dCf_dD(2:end  ).*sT2)                        + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dD2;    
dCDX_dT2= upw    .*(dCD_dT(2:end  ).*sT2 - CD(2:end  ).*sT2./T2) + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dT2;
dCDX_dD2= upw    .*(dCD_dD(2:end  ).*sT2)                        + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dD2; 
dCFX_dU2= upw    .*(dCf_dU(2:end  ).*sT2);
dCDX_dU2= upw    .*(dCD_dU(2:end  ).*sT2);
dCDX_dCt2=upw    .*dCD_dCt(2:end  ).*sT2;

dCFX_dT1= (I-upw).*(dCf_dT(1:end-1).*sT1 - Cf(1:end-1).*sT1./T1) + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dT1;
dCFX_dD1= (I-upw).*(dCf_dD(1:end-1).*sT1)                        + ( Cf(2:end).*sT2 - Cf(1:end-1).*sT1 ) .*dupw_dD1;    
dCDX_dT1= (I-upw).*(dCD_dT(1:end-1).*sT1 - CD(1:end-1).*sT1./T1) + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dT1;
dCDX_dD1= (I-upw).*(dCD_dD(1:end-1).*sT1)                        + ( CD(2:end).*sT2 - CD(1:end-1).*sT1 ) .*dupw_dD1; 
dCFX_dU1= (I-upw).*(dCf_dU(1:end-1).*sT1);
dCDX_dU1= (I-upw).*(dCD_dU(1:end-1).*sT1);
dCDX_dCt1=(I-upw).*dCD_dCt(1:end-1).*sT1;




% shape parameter equation -> right hand side of Newton System
f2= HSlog + (1-HM-HWM).*Ulog + slog.*(CFX - CDX) + h.*VbM./UM.*(1-1./HSM) ./TM;


% derivates of shape parameter equation
%-------------------------------------------------------------------------

% momentum thickness T
df2_dT2 =  1./HS(2:end  ).*dHS_dT(2:end  ) - (0.5*dH_dT(2:end  ) + dHWM_dT(2:end  )).*Ulog ...
            + slog.*(dCFX_dT2 - dCDX_dT2 ) - 0.5*h.*VbM./UM.*(1-1./HSM)./TM.^2 + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dT2;
 
df2_dT1 = -1./HS(1:end-1).*dHS_dT(1:end-1) - (0.5*dH_dT(1:end-1) + dHWM_dT(1:end-1)).*Ulog ...
            + slog.*(dCFX_dT1 - dCDX_dT1 ) - 0.5*h.*VbM./UM.*(1-1./HSM)./TM.^2 + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dT1;
 
 
% displacement thickness DI
df2_dD2 = 1./HS(2:end  ).*dHS_dD(2:end  ) - 0.5*dH_dD(2:end  ).*Ulog + slog.*(dCFX_dD2 - dCDX_dD2) + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dD2 ;

df2_dD1 =-1./HS(1:end-1).*dHS_dD(1:end-1) - 0.5*dH_dD(1:end-1).*Ulog + slog.*(dCFX_dD1 - dCDX_dD1) + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dD1 ;

% shear stress c
df2_dCt2= - slog.*dCDX_dCt2 ;    
df2_dCt1= - slog.*dCDX_dCt1 ; 

% Boundary layer edge velocity U
df2_dU2 =  1./HS(2:end  ).*dHS_dU(2:end  ) + (1-HM-HWM)./U(2:end  ) ...
             + slog.*(dCFX_dU2 - dCDX_dU2) - 0.5*h.*VbM.*(1-1./HSM)./(UM.^2.*TM ) + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dU2 ;            
df2_dU1 = -1./HS(1:end-1).*dHS_dU(1:end-1) - (1-HM-HWM)./U(1:end-1) ...
             + slog.*(dCFX_dU1 - dCDX_dU1) - 0.5*h.*VbM.*(1-1./HSM)./(UM.^2.*TM ) + h.*VbM./(UM.*TM.*HSM.^2).*dHSM_dU1 ; 

% arc length s
df2_ds2= 1./s2.*(CFX - CDX) + slog.*(dCFX_ds2 - dCDX_ds2) + VbM./UM.*(1-1./HSM) ./TM;
df2_ds1=-1./s1.*(CFX - CDX) + slog.*(dCFX_ds1 - dCDX_ds1) - VbM./UM.*(1-1./HSM) ./TM;
            
            
            
% set up quantities for shear stress equation 

Ctlog=log( Ct(2:end)./Ct(1:end-1) );
Ct2=Ct(2:end); Ct1=Ct(1:end-1);


% upwinding
CtM=(I-upw).*Ct(1:end-1) + upw.*Ct(2:end);
CEM=(I-upw).*CE(1:end-1) + upw.*CE(2:end);
CfM=(I-upw).*Cf(1:end-1) + upw.*Cf(2:end);
HM =(I-upw).*H(1:end-1)  + upw.*H(2:end) ;

dHM_dT2 = upw    .*dH_dT(2:end  )  + (H(2:end)-H(1:end-1)).*dupw_dT2;
dHM_dD2 = upw    .*dH_dD(2:end  )  + (H(2:end)-H(1:end-1)).*dupw_dD2;

dCEM_dT2= upw    .*dCE_dT(2:end  ) + (CE(2:end)-CE(1:end-1)).*dupw_dT2;
dCEM_dD2= upw    .*dCE_dD(2:end  ) + (CE(2:end)-CE(1:end-1)).*dupw_dD2;
dCEM_dU2= upw    .*dCE_dU(2:end  );

dCfM_dT2= upw    .*dCf_dT(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT2;
dCfM_dD2= upw    .*dCf_dD(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD2;
dCfM_dU2= upw    .*dCf_dU(2:end  );

dCtM_dT2=(Ct(2:end)-Ct(1:end-1)).*dupw_dT2;
dCtM_dD2=(Ct(2:end)-Ct(1:end-1)).*dupw_dD2;
dCtM_dC2= upw  ;

dHM_dT1 = (I-upw).*dH_dT(1:end-1)  + (H(2:end)-H(1:end-1)).*dupw_dT1;
dHM_dD1 = (I-upw).*dH_dD(1:end-1)  + (H(2:end)-H(1:end-1)).*dupw_dD1;

dCEM_dT1= (I-upw).*dCE_dT(1:end-1) + (CE(2:end)-CE(1:end-1)).*dupw_dT1;
dCEM_dD1= (I-upw).*dCE_dD(1:end-1) + (CE(2:end)-CE(1:end-1)).*dupw_dD1;
dCEM_dU1= (I-upw).*dCE_dU(1:end-1);

dCfM_dT1= (I-upw).*dCf_dT(1:end-1) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT1;
dCfM_dD1= (I-upw).*dCf_dD(1:end-1) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD1;
dCfM_dU1= (I-upw).*dCf_dU(1:end-1);

dCtM_dT1=(Ct(2:end)-Ct(1:end-1)).*dupw_dT1;
dCtM_dD1=(Ct(2:end)-Ct(1:end-1)).*dupw_dD1;
dCtM_dC1= (I-upw);


 % delta with Green correlation
Del=T.*( 3.15 +1.72./(H-1) ) + DI;
dDel_dH = -1.72*T./(H-1).^2;
dDel_dT = (3.15 + 1.72./(H-1)) + dDel_dH.*dH_dT;
dDel_dD = ones(size(DI))       + dDel_dH.*dH_dD;   
    
ind=find(Del > 12*T); % Filter big values of Del
Del(ind)=12*T(ind);  
dDel_dT(ind)=12;
dDel_dD(ind)=0;


DelM=0.5*( Del(1:end-1) + Del(2:end) );


% Uq term in shear stress Equation UQ = 4/3 ( cf/2 - (H-1-18/Re)^2/(6.7*H)^2 )/D
%----------------------------------------------------------
if IsWake
    HKK= HM - 1;
    dHKK_dHM  = ones(dimMi);
    dHKK_dRetM= zeros(dimMi);
    k=0.9; % increase dissipation length on wake
else
    HKK= HM - 1 - 18./RetM;
    dHKK_dHM  = ones(dimMi);
    dHKK_dRetM= 18./RetM.^2;
    ind=find(HKK<0.01);
    HKK(ind)=0.01;
    dHKK_dHM(ind)  =0;
    dHKK_dRetM(ind)=0;
    k=1;
end

HR= HKK./(6.7*k*HM);
dHR_dHKM = dHKK_dHM  ./(6.7*k*HM) - HR./HM;
dHR_dRetM= dHKK_dRetM./(6.7*k*HM);


UQ=(CfM - HR.^2)./(0.75*DM);
dUQ_dHM  =-2*HR.*dHR_dHKM  ./(0.75*DM);
dUQ_dRetM=-2*HR.*dHR_dRetM ./(0.75*DM);
dUQ_dCfM = 1./(0.75*DM);
dUQ_dDM  =-UQ./DM;

dUQ_dT2= dUQ_dHM.*dHM_dT2 + dUQ_dCfM.*dCfM_dT2 ;%+ 0.5*dUQ_dRetM.*dRet_dT(2:end  ) ;
dUQ_dD2= dUQ_dHM.*dHM_dD2 + dUQ_dCfM.*dCfM_dD2 + 0.5*dUQ_dDM                     ;
dUQ_dU2=                    dUQ_dCfM.*dCfM_dU2 ;%+ 0.5*dUQ_dRetM.*dRet_dU(2:end  ) ;

dUQ_dT1= dUQ_dHM.*dHM_dT1 + dUQ_dCfM.*dCfM_dT1 ;%+ 0.5*dUQ_dRetM.*dRet_dT(1:end-1) ;
dUQ_dD1= dUQ_dHM.*dHM_dD1 + dUQ_dCfM.*dCfM_dD1 + 0.5*dUQ_dDM                     ;
dUQ_dU1=                    dUQ_dCfM.*dCfM_dU1 ;%+ 0.5*dUQ_dRetM.*dRet_dU(1:end-1) ;


% factor of the (CtauEQ - Ctau) term 
%----------------------------------------------------------
SC=7.4648./(1+USM);
dSC_dUSM= -SC./(1+USM);
dSC_dT2= 0.5* dSC_dUSM.*dUs_dT(2:end);
dSC_dD2= 0.5* dSC_dUSM.*dUs_dD(2:end);
dSC_dU2= 0.5* dSC_dUSM.*dUs_dU(2:end);
dSC_dT1= 0.5* dSC_dUSM.*dUs_dT(1:end-1);
dSC_dD1= 0.5* dSC_dUSM.*dUs_dD(1:end-1);
dSC_dU1= 0.5* dSC_dUSM.*dUs_dU(1:end-1);

% shear stress Equation ->  Ct=sqrt(Ctau) root already incorporated

f3=SC.*h.*(CEM - k*CtM) + 2*DelM.*(h.*UQ -Ulog -Ctlog) ;
    
% derivates of shear stress equation
%-------------------------------------------------------------------------

% momentum thickness T
df3_dT2 = dSC_dT2.*h.*(CEM - k*CtM) + SC.*h.*(dCEM_dT2 - k*dCtM_dT2) + dDel_dT(2:end  ).*(h.*UQ -Ulog -Ctlog)...
            + 2*DelM.*h.*dUQ_dT2 ;
        
df3_dT1 = dSC_dT1.*h.*(CEM - k*CtM) + SC.*h.*(dCEM_dT1 - k*dCtM_dT1) + dDel_dT(1:end-1).*(h.*UQ -Ulog -Ctlog)...
            + 2*DelM.*h.*dUQ_dT1 ;    
        
% displacement thickness DI
df3_dD2 = dSC_dD2.*h.*(CEM - k*CtM) + SC.*h.*(dCEM_dD2 - k*dCtM_dD2) + dDel_dD(2:end  ).*(h.*UQ -Ulog -Ctlog)...
            + 2*DelM.*h.*dUQ_dD2 ;
        
df3_dD1 = dSC_dD1.*h.*(CEM - k*CtM) + SC.*h.*(dCEM_dD1 - k*dCtM_dD1) + dDel_dD(1:end-1).*(h.*UQ -Ulog -Ctlog)...
            + 2*DelM.*h.*dUQ_dD1 ;             
        
% shear stress c
df3_dCt2= -k*SC.*h.*dCtM_dC2 - 2*DelM./Ct2 ;  
df3_dCt1= -k*SC.*h.*dCtM_dC1 + 2*DelM./Ct1 ; 

% Boundary layer edge velocity U

df3_dU2 = dSC_dU2.*h.*(CEM - k*CtM) + SC.*h.*dCEM_dU2  + 2*DelM.*( h.*dUQ_dU2 -1./U(2:end  ) );
df3_dU1 = dSC_dU1.*h.*(CEM - k*CtM) + SC.*h.*dCEM_dU1  + 2*DelM.*( h.*dUQ_dU1 +1./U(1:end-1) ); 
% arc length s
df3_ds2= SC.*(CEM - k*CtM) + 2*DelM.*UQ;
%df3_ds1=-SC.*(CEM - k*CtM) - 2*DelM.*UQ;
df3_ds1=-df3_ds2;

% function output 

der.df1_dT=[df1_dT1,df1_dT2];
der.df1_dD=[df1_dD1,df1_dD2];
der.df1_dU=[df1_dU1,df1_dU2];
der.df1_dC=[0,0];

der.df2_dT=[df2_dT1,df2_dT2];
der.df2_dD=[df2_dD1,df2_dD2];
der.df2_dU=[df2_dU1,df2_dU2];
der.df2_dC=[df2_dCt1,df2_dCt2];

der.df3_dT=[df3_dT1,df3_dT2];
der.df3_dD=[df3_dD1,df3_dD2];
der.df3_dU=[df3_dU1,df3_dU2];
der.df3_dC=[df3_dCt1,df3_dCt2];

der.df1_ds=[df1_ds1,df1_ds2];
der.df2_ds=[df2_ds1,df2_ds2];
der.df3_ds=[df3_ds1,df3_ds2];
end

