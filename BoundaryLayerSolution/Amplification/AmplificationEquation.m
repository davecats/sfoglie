function [ f, df_dn,df_dT,df_dD,df_dU, df_ds ] = AmplificationEquation(flo,n,T,U,H,Ret,h  )
%AMPLIFICATIONEQUATION evaluates the amplification equation and calculate
%its derivates.

% Calculate the approximated dn/ds
[dn,ddn_dT,ddn_dH, ddn_dRet,ddn_dn] = AmplificationDerivate(flo,H,Ret,T, true, n );

% Amplification equation
% (n2-n1)/h = (dn/ds)_approx
f= n(2:end) - n(1:end-1) - dn.*h;
 
% derivates in respect to amplification exponent n
df_dn1=-ones(size(n(2:end))) - h.*ddn_dn ;
df_dn2= ones(size(n(2:end))) - h.*ddn_dn ;

df_dn=[df_dn1,df_dn2];

% other derivates
dH_dT= -H./T;
dH_dD= 1./T;
dRet_dT= Ret./T;
dRet_dU= Ret./U;

df_dT1= -h.* ( ddn_dT(:,1) + ddn_dH(:,1).*dH_dT(1:end-1) + ddn_dRet(:,1).*dRet_dT(1:end-1) ); 
df_dT2= -h.* ( ddn_dT(:,2) + ddn_dH(:,2).*dH_dT(2:end  ) + ddn_dRet(:,2).*dRet_dT(2:end  ) ); 

df_dD1= -h.*ddn_dH(:,1).*dH_dD(1:end-1); 
df_dD2= -h.*ddn_dH(:,2).*dH_dD(2:end  ); 

df_dU1= -h.*ddn_dRet(:,1).*dRet_dU(1:end-1); 
df_dU2= -h.*ddn_dRet(:,2).*dRet_dU(2:end  ); 

df_dT=[df_dT1,df_dT2];
df_dD=[df_dD1,df_dD2];
df_dU=[df_dU1,df_dU2];
df_ds=[dn,-dn];

end

