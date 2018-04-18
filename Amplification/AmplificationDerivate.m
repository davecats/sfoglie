function [dn,ddn_dT,ddn_dH, ddn_dRet,ddn_dn] = AmplificationDerivate(flo,H,Ret,T,getDerivates,n )
%AMPLIFICATIONEQ    calculates the Amplification increase of each interval dn
%                   and its partial derivates in respect to T, H12 and Ret 

nkrit=flo.nkrit;

if nargin==5
    getDerivates=false;
end

dimP=size(H);
dim=size(H(1:end-1));

k=1./(H-1);
B= 2.492*(k).^0.43;
C=tanh(14*k-9.24);

Gkrit= B + 0.7*(C+1);
GR=log10(Ret);


% zero Amplifikation vor GR < Gkrit
ind=find(GR > Gkrit-0.08 );

RN=(GR(ind)-(Gkrit(ind)-0.08))/0.16;

%----------------- 
Rf=ones(size(RN));
ind2=find(RN<1);
Rtmp=RN(ind2);
Rf(ind2)=3*Rtmp.^2 - 2*Rtmp.^3;
%-----------------    

Arg=3.87.*k(ind)-2.52;
EX= exp(-Arg.^2);
DA=0.028*(H(ind)-1) - 0.0345*EX;
%BRG=-20*k(ind);
Af=-0.05 + 2.7*k(ind)    - 5.5*k(ind).^2 + 3*k(ind).^3 ;%+ 0.1*exp(BRG);


A= zeros(dimP);
A(ind)= (Af.*DA./T(ind)).*Rf;

% add up the amplification delta
Am=sqrt( 0.5*(A(1:end-1).^2 + A(2:end).^2) ); % midpoint approximation

% addition to make shure that dn/ds>0
if getDerivates==false || nargin==5
    n=A;
    %n=zeros(size(A));
end

ARG2=min(20*(nkrit - 0.5*(n(1:end-1)+n(2:end) )) , 20);
EXN=ones(dim);

ind3=find(ARG2>0);
EXN(ind3)=exp(-ARG2(ind3));

add= EXN*0.002./(T(1:end-1)+T(2:end));

dn= (Am + add);

if getDerivates  
    dGkrit_dH= -0.43*B.*k -9.8*k.^2.*(1-C.^2) ;
    dGR_dRet= 1./(2.3025851*Ret);
    
    % RN=(GR-Gkrit-0.08)/0.16;
    dRN_dH  =-dGkrit_dH(ind)/0.16;
    dRN_dRet= dGR_dRet(ind) /0.16;
    
    dRf_dH=zeros(size(RN));
    dRf_dRet=zeros(size(RN));
    
    % ind2 -> RN<1
    dRf_dRN =6*Rtmp.*(1-Rtmp);
    dRf_dH(ind2)  = dRf_dRN.*dRN_dH(ind2);
    dRf_dRet(ind2)= dRf_dRN.*dRN_dRet(ind2);
    
    %Arg=3.87.*k(ind)-2.52;
    dArg_dH=-3.87.*k(ind).^2;

    %EX= exp(-Arg.^2);
    dEX_dH=-2*Arg.*EX.*dArg_dH;

    dDA_dH=0.028 - 0.0345*dEX_dH;

    %BRG=-20*k(ind);
    dAf_dH=  - 2.7*k(ind).^2 + 11 *k(ind).^3 - 9*k(ind).^4 ;% - 2*exp(BRG);
    
    
    dA_dT=zeros(dimP);
    dA_dRet=zeros(dimP);
    dA_dH=zeros(dimP);
    %A(ind)= (Af.*DA./T(ind)).*Rf;
    dA_dH(ind)= (dAf_dH.*DA./T(ind) + Af.*dDA_dH./T(ind)).*Rf + dRf_dH.*Af.*DA./T(ind);
    dA_dT(ind)=-A(ind)./T(ind);
    dA_dRet(ind)= (Af.*DA./T(ind)).*dRf_dRet;
    
    dAm_dA=zeros(dim);
    dAm_dA(Am>0)=1./(2*Am(Am>0));
    
    dAm_dA1= A(1:end-1).*dAm_dA;
    dAm_dA2= A(2:end  ).*dAm_dA;
    
    dEXN_dn=zeros(dim);
    dEXN_dn(ind3)=10*EXN(ind3);

    dadd_dn = dEXN_dn.*0.002./(T(1:end-1)+T(2:end));
    dadd_dT = -EXN./(T(1:end-1)+T(2:end));

    %derivates
    ddn_dA1=  dAm_dA1;
    ddn_dA2=  dAm_dA2;
    
    ddn_dT1= ddn_dA1.*dA_dT(1:end-1) + dadd_dT;
    ddn_dT2= ddn_dA2.*dA_dT(2:end  ) + dadd_dT;

    ddn_dH1= ddn_dA1 .*dA_dH(1:end-1) ;
    ddn_dH2= ddn_dA2 .*dA_dH(2:end  ) ;

    ddn_dRet1= ddn_dA1 .*dA_dRet(1:end-1);
    ddn_dRet2= ddn_dA2 .*dA_dRet(2:end  );

    ddn_dn=dadd_dn;
    ddn_dT=[ddn_dT1,ddn_dT2];
    ddn_dH=[ddn_dH1,ddn_dH2];
    ddn_dRet=[ddn_dRet1,ddn_dRet2];

    
    
%     dEXN_dA=zeros(dim);
%     dEXN_dA(ind3)=10*EXN(ind3);
% 
%     dadd_dA = dEXN_dA.*0.002./(T(1:end-1)+T(2:end));
%     dadd_dT = -EXN./(T(1:end-1)+T(2:end));
% 
%     %derivates
%     ddn_dA1=  (dAm_dA1 + dadd_dA);
%     ddn_dA2=  (dAm_dA2 + dadd_dA);
%     
%     ddn_dT1= ddn_dA1.*dA_dT(1:end-1) + dadd_dT;
%     ddn_dT2= ddn_dA2.*dA_dT(2:end  ) + dadd_dT;
% 
%     ddn_dH1= ddn_dA1 .*dA_dH(1:end-1) ;
%     ddn_dH2= ddn_dA2 .*dA_dH(2:end  ) ;
% 
%     ddn_dRet1= ddn_dA1 .*dA_dRet(1:end-1);
%     ddn_dRet2= ddn_dA2 .*dA_dRet(2:end  );
% 
%     ddn_dA=[ddn_dA1,ddn_dA2];
%     ddn_dT=[ddn_dT1,ddn_dT2];
%     ddn_dH=[ddn_dH1,ddn_dH2];
%     ddn_dRet=[ddn_dRet1,ddn_dRet2];
else
    ddn_dn=zeros(dim);
    ddn_dT=zeros(dim);
    ddn_dH=zeros(dim);
    ddn_dRet=zeros(dim);
end


end

