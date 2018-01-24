function [dn ,dn_dT,dn_dH,dn_dRet] = AmplificationDerivate(H,Ret,T )
%AMPLIFICATIONEQ    calculates the Amplification increase of each interval dn
%                   and its partial derivates in respect to T, H12 and Ret 

dimP=size(H);
dim=size(H(1:end-1));

k=1./(H-1);

B= 2.492*(k).^0.43;


C=tanh(14*k-9.24);

Gkrit= B + 0.7*(C+1);
dGkrit_dH= -0.43*B.*k -9.8*k.^2.*(1-C.^2) ;


GR=log10(Ret);
dGR_dRet= 1./(2.3025851*Ret);


ind=find(GR > Gkrit-0.08 );

A= zeros(dimP);
dA_dH= zeros(dimP);
dA_dT=zeros(dimP);
dA_dRet= zeros(dimP);


if isempty(ind)==false
    RN=(GR(ind)-(Gkrit(ind)-0.08))/0.16;
    dRN_dH  =-dGkrit_dH(ind)/0.16;
    dRN_dRet= dGR_dRet(ind)/0.16;

    Rf=ones(size(RN));
    dRf_dH=zeros(size(RN));
    dRf_dRet=zeros(size(RN));

    ind2=find(RN<1);
    if isempty(ind2)==false
        Rtmp=RN(ind2);
        Rf(ind2)=3*Rtmp.^2 - 2*Rtmp.^3;
        dRf_dRN =6*Rtmp.*(1-Rtmp);

        dRf_dH(ind2)  = dRf_dRN.*dRN_dH(ind2);
        dRf_dRet(ind2)= dRf_dRN.*dRN_dRet(ind2);
    end

    Arg=3.87.*k(ind)-2.52;
    dArg_dH=-3.87.*k(ind).^2;

    EX= exp(-Arg.^2);
    dEX_dH=-2*Arg.*EX.*dArg_dH;

    DA=0.028*(H(ind)-1) - 0.0345*EX;
    dDA_dH=0.028 - 0.0345*dEX_dH;

    %BRG=-20*k(ind);
    Af=-0.05 + 2.7*k(ind)    - 5.5*k(ind).^2 + 3*k(ind).^3 ;%+ 0.1*exp(BRG);
    dAf_dH=  - 2.7*k(ind).^2 + 11 *k(ind).^3 - 9*k(ind).^4 ;% - 2*exp(BRG);

    A(ind)= (Af.*DA./T(ind)).*Rf;
    dA_dH(ind)= (dAf_dH.*DA./T(ind) + Af.*dDA_dH./T(ind)).*Rf + dRf_dH.*Af.*DA./T(ind);
    dA_dT(ind)=-A(ind)./T(ind);
    dA_dRet(ind)= (Af.*DA./T(ind)).*dRf_dRet;
end

% add up the amplification delta

Am=sqrt( 0.5*(A(1:end-1).^2 + A(2:end).^2) ); % midpoint approximation

% prevent singularity
dAm_darg=zeros(dim);
dAm_darg(Am>0)=1./(2*Am(Am>0));

dAm_dA1= dAm_darg .*A(1:end-1);
dAm_dA2= dAm_darg .*A(2:end  );


ARG=min(20*(9-0.5*(A(1:end-1)+A(2:end) )) , 20);

EXN=ones(dim);
dEXN_dA=zeros(dim);

ind=find(ARG>0);

EXN(ind)=exp(-ARG);
dEXN_dA(ind)=10*EXN(ind);


add= EXN*0.002./(T(1:end-1)+T(2:end));
dadd_dA=dEXN_dA*0.002./(T(1:end-1)+T(2:end));
dadd_dH1= dadd_dA.*dA_dH(1:end-1);
dadd_dH2= dadd_dA.*dA_dH(2:end  );
dadd_dRet1= dadd_dA.*dA_dRet(1:end-1);
dadd_dRet2= dadd_dA.*dA_dRet(2:end  );

dadd_dT= -add./(T(1:end-1)+T(2:end));
dadd_dT1= dadd_dT + dA_dT(1:end-1);
dadd_dT2= dadd_dT + dA_dT(2:end  );

dn= (Am+add);

%derivates
dn_dT1=  dAm_dA1.*dA_dT(1:end-1) + dadd_dT1;
dn_dT2=  dAm_dA2.*dA_dT(2:end  ) + dadd_dT2;

dn_dH1= dAm_dA1 .*dA_dH(1:end-1) + dadd_dH1;
dn_dH2= dAm_dA2 .*dA_dH(2:end  ) + dadd_dH2;

dn_dRet1= dAm_dA1 .*dA_dRet(1:end-1) + dadd_dRet1;
dn_dRet2= dAm_dA2 .*dA_dRet(2:end  ) + dadd_dRet2;


dn_dT=[dn_dT1,dn_dT2];
dn_dH=[dn_dH1,dn_dH2];
dn_dRet=[dn_dRet1,dn_dRet2];

end

