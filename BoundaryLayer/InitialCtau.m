function [C,dC_dT,dC_dD,dC_dU]= InitialCtau( D,T,U,nu, getDer )
%INITIALCTAU  returns the initial values for the maximum shear stress
%             coefficient Ctau at the first turbulent node



H=D./T;
dH_dT=-H./T;
dH_dD=1./T;

Ret=U.*T/nu;
dRet_dT = U/nu;
dRet_dU = T/nu;

% empirical factor 
fac=1.8*exp(-3.3./(H-1));
%fac=1.1*exp(-10./HK2.^2);

if nargin==4 || getDer==false
    
    [ HS, ~,~   ]=H32turb( H,Ret);
    Us=0.5*HS.*( -1/3 + T./(0.75*D) );

    % equilibrium Ctau
    %
    [CEQ,~,~,~, ~]=CtEQ( H,Ret,HS,Us,false);
    C=fac.* CEQ;

else
    [ HS, dHS_dH,dHS_dRet   ]=H32turb( H,Ret);
    Us=0.5*HS.*( -1/3 + T./(0.75*D) );
    
    [CEQ,dCE_dH,dCE_dRet,dCE_dHS, dCE_dUs]=CtEQ( H,Ret,HS,Us,false);
    
    
    C=fac.* CEQ;
    
    dfac_dH= 3.3./(H-1).^2 .*fac;
    dfac_dD=dfac_dH.*dH_dD;
    dfac_dT=dfac_dH.*dH_dT;


    dHS_dT  = dHS_dH.*dH_dT + dHS_dRet.*dRet_dT;
    dHS_dD  = dHS_dH.*dH_dD;
    dHS_dU  =                dHS_dRet.*dRet_dU;

    dUs_dHS=  Us./HS;
    dUs_dH = -HS./(1.5*H.^2);

    dUs_dT = dUs_dH.*dH_dT + dUs_dHS.*dHS_dT;
    dUs_dD = dUs_dH.*dH_dD + dUs_dHS.*dHS_dD;
    dUs_dU =                 dUs_dHS.*dHS_dU;


    dCE_dT = dCE_dH.*dH_dT + dCE_dRet.*dRet_dT + dCE_dHS.*dHS_dT + dCE_dUs.*dUs_dT ;
    dCE_dD = dCE_dH.*dH_dD +                   + dCE_dHS.*dHS_dD + dCE_dUs.*dUs_dD ;
    dCE_dU =                 dCE_dRet.*dRet_dU + dCE_dHS.*dHS_dU + dCE_dUs.*dUs_dU ;

    dC_dT=dfac_dT.* CEQ + fac.*dCE_dT;
    dC_dD=dfac_dD.* CEQ + fac.*dCE_dD;
    dC_dU= fac.*dCE_dU;
end

end

