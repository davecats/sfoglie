function [sT, der ] = FindTransitionpoint( flo, n, T, D,U,s1,h, mode )
%FINDTRANSITIONPOINT Finds the exakt transition point on the panel where
%transition occurs (n2>nkrit && n1<nkrit) and calculate its derivates in
%respect to the variables U,D,T and n
% mode: defines the approximative differnce scheme for the derivate dn/ds used on the transition panel
%       mode=1: (n2-n1)/(s2-s1) \approx  n'(sT)   -> second order scheme (default)
%       mode=2: (nkrit-n1)/(sT-s1) \approx n'(s1) -> first order scheme (more stable)


nkrit=flo.nkrit;
nu=flo.nu;

H=D./T;
Ret=T.*U/nu;

dH_dD=1./T;
dH_dT=-H./T;

dRet_dT=U/nu;
dRet_dU=T/nu;

s2=s1+h;

if nargin==7
   mode=1; 
end


if mode == 1
    % second order method of Xfoil
    %   equation 
    %     (n2-n1)/(s2-s1) \approx  n'(sT)
    %   is solved for sT. 
    
    k=0; dn_it=1;
    
    
    while abs(dn_it)>5e-5 && k<30
        % weigthening factors between 1 and 2 node   
        
        w2= (nkrit-n(1))/(n(2)-n(1));
        dw2_dn2= -w2/(n(2)-n(1));
        dw2_dn1= (w2-1)/(n(2)-n(1));
        
        w1= 1-w2;
        dw1_dn2= -dw2_dn2;
        dw1_dn1= -dw2_dn1;
        
        % Transition point quantities
        TT=w1*T(1) + w2*T(2);
        UT=w1*U(1) + w2*U(2);
        DT=w1*D(1) + w2*D(2);
        sT=w1*s1   + w2*s2;
        
        % derivates in respect to n2
        dTT_dn2= T(1)*dw1_dn2 + T(2)*dw2_dn2;
        dDT_dn2= D(1)*dw1_dn2 + D(2)*dw2_dn2;
        dUT_dn2= U(1)*dw1_dn2 + U(2)*dw2_dn2;
        dsT_dn2= s1  *dw1_dn2 + s2  *dw2_dn2 ;

        HT=max(1.05,DT/TT);
        RetT= UT*TT/nu;
        
        dHT_dD=1/TT;
        dHT_dT=-HT/TT;

        dRetT_dT=UT/nu;
        dRetT_dU=TT/nu;
        
        % n'(sT)
        [dn,ddn_dT,ddn_dH, ddn_dRet,ddn_dn ] = AmplificationDerivate(flo,[H(1); HT],[Ret(1); RetT],[T(1);TT],true,[n(1) ;nkrit]);
        
        if dn<0; break; end % leave loop in case there is a Amplification decrease
        %dn=max(dn,0);
        
        % derivate of n'(sT) in resprect to 
        ddn_dn2=   ( ddn_dT(2) + ddn_dH(2)*dHT_dT + ddn_dRet(2)*dRetT_dT)*dTT_dn2 ...
                 + (             ddn_dH(2)*dHT_dD                       )*dDT_dn2 ... 
                 + (                                ddn_dRet(2)*dRetT_dU)*dUT_dn2 ;
        
         %(n2-n1)/(s2-s1) \approx  n'(sT) = g(w1,w2, T1,T2,...)  
         % find n2 that satifys the condition
        f=n(2)-n(1)-h*dn;
        
        df_dn2=1-h*ddn_dn2;

        dn_it= -f./df_dn2;
        DsT=dn_it*dsT_dn2;

        rel=1; % relaxation
        if abs(DsT/h)>0.05; rel=0.05*abs(h/DsT);end
        if rel*abs(dn_it)>1; rel=1/abs(dn_it);end

        nT=n(2) + rel*dn_it;
        
        if (n(2)>nkrit && nT<nkrit) || (n(2)<nkrit && nT>nkrit)
            nT=nkrit;
        end

        n(2)=nT;
        
        k=k+1;
    end
    
    % if weigthening factors don�t make sense -> try 1. order method
    if w1>1 || w2>1 || w1<0 ||w2<0
        [sT, der ] = FindTransitionpoint( flo, n, T, D,U,s1,h, 2 );
        return;
    end
    
    
    % derivates of the Transition point quantities
    dsT_dn1= s1  * dw1_dn1 + s2  * dw2_dn1;
    dTT_dn1= T(1)* dw1_dn1 + T(2)* dw2_dn1;
    dDT_dn1= D(1)* dw1_dn1 + D(2)* dw2_dn1;
    dUT_dn1= U(1)* dw1_dn1 + U(2)* dw2_dn1;
    dsT_dn2= s1  * dw1_dn2 + s2  * dw2_dn2;
    dTT_dn2= T(1)* dw1_dn2 + T(2)* dw2_dn2;
    dDT_dn2= D(1)* dw1_dn2 + D(2)* dw2_dn2;
    dUT_dn2= U(1)* dw1_dn2 + U(2)* dw2_dn2;


    % derivates of the Amplification increase dn= nT - n1
    ddn_dT1=     ddn_dT(1) + ddn_dH(1)*dH_dT(1) + ddn_dRet(1)*dRet_dT(1) ...
            + ( ddn_dT(2) + ddn_dH(2)*dHT_dT   + ddn_dRet(2)*dRetT_dT )*w1;

    ddn_dD1=     ddn_dH(1)*dH_dD(1) +  ddn_dH(2)*dHT_dD *w1;
    ddn_dU1=     ddn_dRet(1)*dRet_dU(1) + ddn_dRet(2)*dRetT_dT *w1;

    ddn_dn1=   ( ddn_dT(2) + ddn_dH(2)*dHT_dT + ddn_dRet(2)*dRetT_dT)*dTT_dn1 ...
             + (             ddn_dH(2)*dHT_dD                       )*dDT_dn1 ... 
             + (                                ddn_dRet(2)*dRetT_dU)*dUT_dn1 ...
             + ddn_dn   ;


    ddn_dT2=   ( ddn_dT(2) + ddn_dH(2)*dHT_dT   + ddn_dRet(2)*dRetT_dT )*w2;
    ddn_dD2=    ddn_dH(2)*dHT_dD *w2;
    ddn_dU2=    ddn_dRet(2)*dRetT_dT *w2;



    ddn_dn2=   ( ddn_dT(2) + ddn_dH(2)*dHT_dT + ddn_dRet(2)*dRetT_dT)*dTT_dn2 ...
             + (             ddn_dH(2)*dHT_dD                       )*dDT_dn2 ... 
             + (                                ddn_dRet(2)*dRetT_dU)*dUT_dn2 ;

    % derivates of the Amplification equation
    df_dn1= -1 - h*ddn_dn1;
    df_dT1= -h*ddn_dT1;
    df_dD1= -h*ddn_dD1;
    df_dU1= -h*ddn_dU1;
    df_ds1=  dn;

    df_dn2=  1 - h*ddn_dn2;
    df_dT2= -h*ddn_dT2;
    df_dD2= -h*ddn_dD2;
    df_dU2= -h*ddn_dU2;
    df_ds2= -dn;

    % derivates of the transition arclength sT depending on 1 and 2 values
    der.dsT_dn1 = dsT_dn1 - dsT_dn2/ df_dn2 * df_dn1;
    der.dsT_dT1 =         - dsT_dn2/ df_dn2 * df_dT1;
    der.dsT_dD1 =         - dsT_dn2/ df_dn2 * df_dD1;
    der.dsT_dU1 =         - dsT_dn2/ df_dn2 * df_dU1;
    der.dsT_ds1 =  w1     -dsT_dn2/ df_dn2 * df_ds1;

    der.dsT_dT2 =         - dsT_dn2/ df_dn2 * df_dT2;
    der.dsT_dD2 =         - dsT_dn2/ df_dn2 * df_dD2;
    der.dsT_dU2 =         - dsT_dn2/ df_dn2 * df_dU2;
    der.dsT_ds2 =  w2     -dsT_dn2/ df_dn2 * df_ds2;
    
%--------------------------------------------------------------------------------------------------    
    
elseif mode==2    
    % first order method of Xfoil
    %   ->  (nkrit-n1)/(sT-s1) \approx n'(s1)
    %    
    
    % n'(sT)
    [dn,ddn_dT,ddn_dH, ddn_dRet,ddn_dn] = AmplificationDerivate(flo,H,Ret,T,true,n);
    
    sT=s1 + (nkrit-n(1))/dn;
    
    dsT_ddn=(n(1)-nkrit)/dn^2;

    % derivates of sT in respect to variables at node 1 and 2
    der.dsT_dn1 = -1/dn + dsT_ddn*ddn_dn(1) ;
    der.dsT_ds1 = 1;
    der.dsT_dT1 = dsT_ddn*( ddn_dT(1) + ddn_dH(1)*dH_dT(1) + ddn_dRet(1)*dRet_dT(1));
    der.dsT_dD1 = dsT_ddn*(             ddn_dH(1)*dH_dD(1)                         );
    der.dsT_dU1 = dsT_ddn*(                                  ddn_dRet(1)*dRet_dU(1));
    
    der.dsT_ds2 = 0;
    der.dsT_dT2 = dsT_ddn*( ddn_dT(2) + ddn_dH(2)*dH_dT(2) + ddn_dRet(2)*dRet_dT(2));
    der.dsT_dD2 = dsT_ddn*(             ddn_dH(2)*dH_dD(2)                         );
    der.dsT_dU2 = dsT_ddn*(                                  ddn_dRet(2)*dRet_dU(2));
    
end



















end

