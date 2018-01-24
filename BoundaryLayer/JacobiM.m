function [ J, rhs] = JacobiM( T,m,U,DI,D ,h)
%JACOBIM        calculates the Jacobimatrix and right hand side for the Newton system 
%               for solution vector z=[T1,..,TN,m1,..,mN]^T 
%               gamma: circulation distribution
%               D: Coefficients of mass Defekt -> U = Uinv + Dm
%               L: Vector with panel length
%               m: mass defect vector from previous step
%               T: momentum thickness vector from previous step


nu=evalin('base','nu');
Nle=evalin('base','Nle');
NW=evalin('base','NW');     % number of wake nodes
N=length(T); 
NFoil=N-NW;% node number on air foil 

% panel midpoint approximations
TM=(T(1:end-1) + T(2:end)) /2;
DM=(DI(1:end-1)+ DI(2:end))/2;
UM=(U(1:end-1) + U(2:end)) /2;


% help quantities
H=DI./T; %shape parameter H1
%HM=(H(2:end)+H(1:end-1));%HM=DM./TM;
Ret= U.*T/nu;
%ReM=(Ret(2:end)+Ret(1:end-1));%UM.*TM/nu;

% set Hmin
H(H(1:NFoil)<1.05)=1.05;
ind=find(H(NFoil+1)<1.00005);ind=NFoil*ones(length(ind))+ind;
H(ind)=1.00005;
 
% calculate Cf, CD, H32 and their derivates in respect to H and Ret
[ Cf, dCf_dH, dCf_dRet ] = CF2lam( H(1:end-NW), Ret(1:end-NW));
[ HS, dHS_dH ]=H32lam(H);
[ CD, dCD_dH, dCD_dRet ]=CD2lam(H(1:end-NW),Ret(1:end-NW),HS(1:end-NW));
dCD_dHS=CD./HS(1:end-NW);

% no friction on wake
Cf=[Cf; zeros(NW,1)];dCf_dH=[dCf_dH; zeros(NW,1)];dCf_dRet=[dCf_dRet; zeros(NW,1)];

[ CDW, dCDW_dH,dCDW_dRet ] = CD2Lwake( H(NFoil+1:end),Ret(NFoil+1:end) );
 % wake has 2 sides of boundary layer -> double Dissipation 
CD=[CD; 2*CDW];dCD_dH=[dCD_dH; 2*dCDW_dH]; dCD_dRet =[dCD_dRet; 2*dCDW_dRet];
dCD_dHS=[dCD_dHS; zeros(NW,1)];

% mid point values
CDM=(CD(2:end) + CD(1:end-1))/2;
HSM=(HS(2:end) + HS(1:end-1))/2;

HM= (H(1:end-1)+H(2:end))/2;
RetM=(Ret(1:end-1)+Ret(2:end))/2;
[CfM,dCfM_dH,dCfM_dRet]=CF2lam( HM, RetM);

CfM= CfM/2 + (Cf(2:end) + Cf(1:end-1))/4;

% Derivates T, D, U
%-------------------------------------------------------------------------

% derivates in respect to T
dRet_dT = U/nu;
dH_dT   =-H./T;
dHS_dT  = dHS_dH.*dH_dT;

dCf_dT= dCf_dH.*dH_dT + dCf_dRet.*dRet_dT;
dCD_dT= dCD_dH.*dH_dT + dCD_dRet.*dRet_dT + dCD_dHS.*dHS_dT;

dCfM_dT1=dCfM_dH.*dH_dT(1:end-1)/4 + dCfM_dRet.*dRet_dT(1:end-1)/4 + dCf_dT(1:end-1)/4;
dCfM_dT2=dCfM_dH.*dH_dT(2:end  )/4 + dCfM_dRet.*dRet_dT(2:end  )/4 + dCf_dT(2:end  )/4;

%derivates in respect to D
dH_dD   = 1./T;
dCf_dD  = dCf_dH .*dH_dD;
dHS_dD = dHS_dH.*dH_dD;
dCD_dD  = dCD_dH .*dH_dD + dCD_dHS.*dHS_dD;

dCfM_dD1=dCfM_dH.*dH_dD(1:end-1)/4 + dCf_dD(1:end-1)/4;
dCfM_dD2=dCfM_dH.*dH_dD(2:end  )/4 + dCf_dD(2:end  )/4;

%derivates in respect to U
dRet_dU= T/nu;
dCf_dU = dCf_dRet.*dRet_dU;
dCD_dU = dCD_dRet.*dRet_dU;

dCfM_dU1=  dCfM_dRet.*dRet_dU(1:end-1)/4 + dCf_dU(1:end-1)/4;
dCfM_dU2=  dCfM_dRet.*dRet_dU(2:end  )/4 + dCf_dU(2:end  )/4;

%-------------------------------------------------------------------------




% on suction side point "1" and "2" are swapped because the stream goes in
% opposit direction to the indexing
sgn=[-ones(Nle-1,1);ones(N-Nle,1)];


% finite differences of each panel
dT =sgn.*(T(2:end) - T(1:end-1) );
dU =sgn.*(U(2:end) - U(1:end-1) );
dHS=sgn.*(HS(2:end)-HS(1:end-1) );




% momentum equation -> right hand side of Newton System
f1=dT./h + (2*TM + DM).*dU./(h.*UM) - CfM;


% derivates of momentum equation
%-------------------------------------------------------------------------

% momentum thickness T
    % partial derivate
    adT1=-sgn./h + dU./(h.*UM);
    adT2= sgn./h + dU./(h.*UM);
    
    % total derivates
    df1_dT1 =  adT1 -  dCfM_dT1;
    df1_dT2 =  adT2 -  dCfM_dT2;
    
    clear a1dT1 a1dT2 
    
% displacement thickness D
    % partial derivate
    adD1= dU./(2*h.*UM) ;
    adD2= dU./(2*h.*UM);
    
    % total derivates
    df1_dD1 =  adD1 -  dCfM_dD1;
    df1_dD2 =  adD2 -  dCfM_dD2;

    clear a1dD1 a1dD2 


% tangential velocity U
    % partial derivate
    adU1=-sgn.*(2*TM+DM)./(h.*UM) -(2*TM+DM).*dU./(2*h.*UM.^2);
    adU2= sgn.*(2*TM+DM)./(h.*UM) -(2*TM+DM).*dU./(2*h.*UM.^2);

    
    % total derivates
    df1_dU1 =  adU1 -  dCfM_dU1;
    df1_dU2 =  adU2 -  dCfM_dU2;

    clear a1dU1 a1dU2  




%-------------------------------------------------------------------------

% upwinding for shape parameter Equation
[upw, dupw_dH1,dupw_dH2]=Upwinding(H(1:NFoil-1),H(2:NFoil),false);
[upwW, dupwW_dH1,dupwW_dH2]=Upwinding(H(NFoil+1:end-1),H(NFoil+2:end),true);
upw=[upw; 0.5;upwW];dupw_dH1=[dupw_dH1;0.5;dupwW_dH1];dupw_dH2=[dupw_dH2;0.5;dupwW_dH2];
dupw_dT1=dupw_dH1.*dH_dT(1:end-1);
dupw_dT2=dupw_dH2.*dH_dT(2:end  );
dupw_dD1=dupw_dH1.*dH_dD(1:end-1);
dupw_dD2=dupw_dH2.*dH_dD(2:end  );


I=ones(N-1,1);
CfM=(I-upw).*Cf(1:end-1) + upw.*Cf(2:end);
CDM=(I-upw).*CD(1:end-1) + upw.*CD(2:end);

dCfM_dT1= (I-upw).*dCf_dT(1:end-1) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT1;
dCfM_dT2= upw    .*dCf_dT(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dT2;
dCfM_dD1= (I-upw).*dCf_dD(1:end-1) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD1;
dCfM_dD2= upw    .*dCf_dD(2:end  ) + (Cf(2:end)-Cf(1:end-1)).*dupw_dD2;
dCDM_dT1= (I-upw).*dCD_dT(1:end-1) + (CD(2:end)-CD(1:end-1)).*dupw_dT1;
dCDM_dT2= upw    .*dCD_dT(2:end  ) + (CD(2:end)-CD(1:end-1)).*dupw_dT2;
dCDM_dD1= (I-upw).*dCD_dD(1:end-1) + (CD(2:end)-CD(1:end-1)).*dupw_dD1;
dCDM_dD2= upw    .*dCD_dD(2:end  ) + (CD(2:end)-CD(1:end-1)).*dupw_dD2;
dCfM_dU1= (I-upw).*dCf_dU(1:end-1) ;
dCfM_dU2= upw    .*dCf_dU(2:end  ) ;
dCDM_dU1= (I-upw).*dCD_dU(1:end-1) ;
dCDM_dU2= upw    .*dCD_dU(2:end  ) ;

% dCfM_dT1= dCf_dT(1:end-1) /2;
% dCfM_dT2= dCf_dT(2:end  ) /2;
% dCfM_dD1= dCf_dD(1:end-1) /2;
% dCfM_dD2= dCf_dD(2:end  ) /2;
% dCDM_dT1= dCD_dT(1:end-1) /2;
% dCDM_dT2= dCD_dT(2:end  ) /2;
% dCDM_dD1= dCD_dD(1:end-1) /2;
% dCDM_dD2= dCD_dD(2:end  ) /2;
% dCfM_dU1= dCf_dU(1:end-1) /2;
% dCfM_dU2= dCf_dU(2:end  ) /2;
% dCDM_dU1= dCD_dU(1:end-1) /2;
% dCDM_dU2= dCD_dU(2:end  ) /2;


% shape parameter equation -> right hand side of Newton System
f2= TM.*dHS./h + HSM.*(TM-DM).*dU./(h.*UM) - CDM + HSM.*CfM;


% derivates of  shape parameter equation 
%-------------------------------------------------------------------------


 % partial derivate in Respect to HS
    df2_dHS1= -sgn.*TM./h + (TM-DM).*dU./(2*h.*UM) + CfM/2;
    df2_dHS2=  sgn.*TM./h + (TM-DM).*dU./(2*h.*UM) + CfM/2;


% momentum thickness T
    % partial derivate
    adT1= dHS./(2*h) + HSM.*dU./(2*h.*UM);
    adT2= dHS./(2*h) + HSM.*dU./(2*h.*UM);

    % total derivates
    df2_dT1 =  adT1 + df2_dHS1.*dHS_dT(1:end-1) - dCDM_dT1 + HSM.*dCfM_dT1;
    df2_dT2 =  adT2 + df2_dHS2.*dHS_dT(2:end)   - dCDM_dT2 + HSM.*dCfM_dT2;
    
    clear a1dT1 a1dT2 
    
% displacement thickness D
    % partial derivate
    adD1= -HSM.*dU./(2*h.*UM) ;
    adD2= -HSM.*dU./(2*h.*UM) ;
    
    % total derivates
    df2_dD1 =  adD1 + df2_dHS1.*dHS_dD(1:end-1) - dCDM_dD1 + HSM.*dCfM_dD1;
    df2_dD2 =  adD2 + df2_dHS2.*dHS_dD(2:end)   - dCDM_dD2 + HSM.*dCfM_dD2;

    clear a1dD1 a1dD2 


% tangential velocity U
    % partial derivate
    adU1= -sgn.*HSM.*(TM-DM)./(h.*UM) - HSM.*(TM-DM).*dU./(2*h.*UM.^2) ;
    adU2=  sgn.*HSM.*(TM-DM)./(h.*UM) - HSM.*(TM-DM).*dU./(2*h.*UM.^2) ;

    
    % total derivates -> dHS_dU=0
    df2_dU1 =  adU1  - dCDM_dU1 + HSM.*dCfM_dU1;
    df2_dU2 =  adU2  - dCDM_dU2 + HSM.*dCfM_dU2;

    clear a1dU1 a1dU2  
    
    
% switch to m as variable

%derivates of D,U in respect to m
dD_dm=diag( 1./U ) - ( m./(U.^2)*ones(1,N) ).* D;
dD1_dm=dD_dm(1:end-1,:);
dD2_dm=dD_dm(2:end,:);

dU1_dm=D(1:end-1,:);
dU2_dm=D(2:end,:);


% total derivates of equations in respect to m   


% equation 1-> momentum Eq
df1_dm=  (df1_dD1*ones(1,N)) .* dD1_dm  +  (df1_dD2*ones(1,N)) .* dD2_dm ...
        +(df1_dU1*ones(1,N)) .* dU1_dm  +  (df1_dU2*ones(1,N)) .* dU2_dm;


% equation 2-> shape parameter Eq
df2_dm=  (df2_dD1*ones(1,N)) .* dD1_dm  +  (df2_dD2*ones(1,N)) .* dD2_dm ...
        +(df2_dU1*ones(1,N)) .* dU1_dm  +  (df2_dU2*ones(1,N)) .* dU2_dm;
    

 %-------------------------------------------------------------------------
 


% create Jacobi matrix for Newton System

% contribution of momentum eq
JT1= [diag(df1_dT1), zeros(N-1,1)]; % last node is never a "1" node
JT2= [zeros(N-1,1), diag(df1_dT2)]; % first node is never a "2" node
JTf1= JT1+JT2;


% contribution of shape parameter eq
JT1= [diag(df2_dT1), zeros(N-1,1)]; % last node is never a "1" node
JT2= [zeros(N-1,1), diag(df2_dT2)]; % first node is never a "2" node
JTf2= JT1+JT2;


%Total T part
JT=zeros(2*length(f1),N);
JT(1:2:end-1,1:1:N)=JTf1;
JT(2:2:end,1:1:N)=JTf2;

Jm=zeros(2*length(f1),N);
Jm(1:2:end-1,:)=df1_dm;
Jm(2:2:end,:)=df2_dm;

J=[JT,Jm];

rhs=zeros(2*length(f1),1);
rhs(1:2:end-1)=f1;
rhs(2:2:end)=f2;

% J(2*Nle-3,:)=[];J(2*Nle-2,:)=[];rhs(2*Nle-3)=[];rhs(2*Nle-2)=[];
% J(2*NFoil-1,:)=[];J(2*NFoil,:)=[];rhs(2*NFoil-1)=[];rhs(2*NFoil)=[];


% Set starting conditions to close LGS
%---------------------------------------------------

% -> values on Leading edge set by initial conditions -> no change
J(2*Nle-3,:)=zeros(1,length(J(1,:))); J(2*Nle-3,Nle-1)=1;
J(2*Nle-2,:)=zeros(1,length(J(1,:))); J(2*Nle-2,Nle)=1;
rhs(2*Nle-3)=0;
rhs(2*Nle-2)=0;
ms=zeros(1,length(J(1,:))); ms(N+Nle-1)=1;
mp=zeros(1,length(J(1,:))); mp(N+Nle)=1;
J=[J; ms;mp];
rhs=[rhs; 0;0];

% values for first wake node set by initial conditions -> no change
J(2*NFoil-1,:)=zeros(1,length(J(1,:))); J(2*NFoil-1,NFoil+1)=1;
J(2*NFoil,:)  =zeros(1,length(J(1,:))); J(2*NFoil,N+NFoil+1)=1;
rhs(2*NFoil-1)=0;
rhs(2*NFoil)  =0;

% minus from Newton method
rhs=-rhs;

end

