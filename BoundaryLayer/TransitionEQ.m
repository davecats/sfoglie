function [ f, erg,sT ] = TransitionEQ( flo,eng, n, T, D,U,Vb,s1,h,C2,sF ,IsForced)
%TRANSITIONEQ   gives the equationsystem for the transition panel
%               -> uses sum of laminar and turbulent part for momentum and shapeparameter equation
%               -> uses Ctau Equation only for turbulent part

nkrit=flo.nkrit;
nu=flo.nu;
TranEQ2=eng.tranEQ;

if nargin==10; IsForced=false;  end

if TranEQ2 && ~IsForced
   [ f, erg,sT ] = TransitionEQ2( flo,eng,n, T, D,U,Vb,s1,h, C2) ; return; 
end


H=D./T;
Ret=T.*U/nu;

dH_dD=1./T;
dH_dT=-H./T;

dRet_dT=U/nu;
dRet_dU=T/nu;

dnI=1;
k=0;

s2=s1+h;

if IsForced
    % set weightening factors in case of forced transition
   w2= (sF-s1)/h ;
   w1= 1-w2;       

   dw2_dn1= 0 ;dw1_dn1= 0 ;
   dw2_dn2= 0 ;dw1_dn2= 0 ;   
   
end


if ~IsForced
    %n2old=n(2);
    
    % solves the Amplification Equation for current step forcing n(sT)= nkrit
    while abs(dnI)>5e-5 && k<30
        w2= (nkrit-n(1))/(n(2)-n(1));
        dw2_dn2= -w2/(n(2)-n(1));
        dw2_dn1= (w2-1)/(n(2)-n(1));
        w1= 1-w2;
        dw1_dn2= -dw2_dn2;
        dw1_dn1= -dw2_dn1;

        TT=w1*T(1) + w2*T(2);
        UT=w1*U(1) + w2*U(2);
        DT=w1*D(1) + w2*D(2);
        sT=w1*s1   + w2*s2;    

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

        [dn,ddn_dT,ddn_dH, ddn_dRet,ddn_dn] = AmplificationDerivate(flo,[H(1); HT],[Ret(1); RetT],[T(1);TT],true,[n(1) ;nkrit]);

        if dn<0; break; end % leave loop in case there is a Amplification decrease
        
        % derivate of amplification increase
        ddn_dn2=   ( ddn_dT(2) + ddn_dH(2)*dHT_dT + ddn_dRet(2)*dRetT_dT)*dTT_dn2 ...
                 + (             ddn_dH(2)*dHT_dD                       )*dDT_dn2 ... 
                 + (                                ddn_dRet(2)*dRetT_dU)*dUT_dn2 ;

        f=n(2)-n(1)-h*dn; % Amplification Equation
        df_dn2=1-h*ddn_dn2;

        dnI= -f./df_dn2;
        DsT=dnI*dsT_dn2;

        rel=1; % relaxation
        if abs(DsT/h)>0.05; rel=0.05*abs(h/DsT);end
        if rel*abs(dnI)>1; rel=1/abs(dnI);end

        nT=n(2) + rel*dnI;
        
        
        if (n(2)>nkrit && nT<nkrit) || (n(2)<nkrit && nT>nkrit)
            nT=nkrit;
        end

        n(2)=nT;
        
        k=k+1;
    end % -> end of iteration
    
end


if IsForced 
    sT= sF;
    % derivates vanish since transition arc length is fix
    dsT_dn1 =0;
    dsT_dT1 =0;
    dsT_dD1 =0;
    dsT_dU1 =0;
    dsT_ds1 =0;
    dsT_dT2 =0;
    dsT_dD2 =0;
    dsT_dU2 =0;
    dsT_ds2 =0;
else
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
    dsT_dn1 = dsT_dn1 - dsT_dn2/ df_dn2 * df_dn1;
    dsT_dT1 =         - dsT_dn2/ df_dn2 * df_dT1;
    dsT_dD1 =         - dsT_dn2/ df_dn2 * df_dD1;
    dsT_dU1 =         - dsT_dn2/ df_dn2 * df_dU1;
    dsT_ds1 =  w1     -dsT_dn2/ df_dn2 * df_ds1;

    dsT_dT2 =         - dsT_dn2/ df_dn2 * df_dT2;
    dsT_dD2 =         - dsT_dn2/ df_dn2 * df_dD2;
    dsT_dU2 =         - dsT_dn2/ df_dn2 * df_dU2;
    dsT_ds2 =  w2     -dsT_dn2/ df_dn2 * df_ds2;
end  
 

%------------------------------------------------------------------------------------------------------
%                       get the boundary layer equations for Transition panel



% Transition point value derivates
   
w2= (sT-s1) /h;
w1=1-w2;

% errormessage in case of wrong weightening factors
% if w1<-1 || w1>2
%    str=['IsForced ',num2str(IsForced) ];
%    if IsForced; str=[str,'sF ', num2str(sF),'s1 ', num2str(s1),'s2 ', num2str(s2) ];end
%    disp(['WARNING:',str]);
%    %error('TransitionEQ weightening factors not in range 0<w<1');
% end

if w2<1e-8; w1=1; w2=0; sT=s1; end
if w1<1e-8; w1=0; w2=1; sT=s2; end





dw2_dsT=1/h;

dw2_dn1= dw2_dsT*dsT_dn1;
dw2_ds1= dw2_dsT*dsT_ds1 + (w2-1) /h;
dw2_ds2= dw2_dsT*dsT_ds2 -  w2    /h;
dw2_dT1= dw2_dsT*dsT_dT1;
dw2_dT2= dw2_dsT*dsT_dT2;
dw2_dD1= dw2_dsT*dsT_dD1;
dw2_dD2= dw2_dsT*dsT_dD2;
dw2_dU1= dw2_dsT*dsT_dU1;
dw2_dU2= dw2_dsT*dsT_dU2;


dw1_dn1= -dw2_dn1;
dw1_ds1= -dw2_ds1;
dw1_ds2= -dw2_ds2;
dw1_dT1= -dw2_dT1;
dw1_dT2= -dw2_dT2;
dw1_dD1= -dw2_dD1;
dw1_dD2= -dw2_dD2;
dw1_dU1= -dw2_dU1;
dw1_dU2= -dw2_dU2;


TT=w1*T(1) + w2*T(2);
dTT_dn1= T(1)*dw1_dn1  + T(2)* dw2_dn1;
dTT_ds1= T(1)*dw1_ds1  + T(2)* dw2_ds1;
dTT_ds2= T(1)*dw1_ds2  + T(2)* dw2_ds2;
dTT_dT1= T(1)*dw1_dT1  + T(2)* dw2_dT1 + w1;
dTT_dT2= T(1)*dw1_dT2  + T(2)* dw2_dT2 + w2;
dTT_dD1= T(1)*dw1_dD1  + T(2)* dw2_dD1;
dTT_dD2= T(1)*dw1_dD2  + T(2)* dw2_dD2;
dTT_dU1= T(1)*dw1_dU1  + T(2)* dw2_dU1;
dTT_dU2= T(1)*dw1_dU2  + T(2)* dw2_dU2;

DT=w1*D(1) + w2*D(2);
dDT_dn1= D(1)*dw1_dn1  + D(2)* dw2_dn1;
dDT_ds1= D(1)*dw1_ds1  + D(2)* dw2_ds1;
dDT_ds2= D(1)*dw1_ds2  + D(2)* dw2_ds2;
dDT_dT1= D(1)*dw1_dT1  + D(2)* dw2_dT1;
dDT_dT2= D(1)*dw1_dT2  + D(2)* dw2_dT2;
dDT_dD1= D(1)*dw1_dD1  + D(2)* dw2_dD1 + w1;
dDT_dD2= D(1)*dw1_dD2  + D(2)* dw2_dD2 + w2;
dDT_dU1= D(1)*dw1_dU1  + D(2)* dw2_dU1;
dDT_dU2= D(1)*dw1_dU2  + D(2)* dw2_dU2;

UT=w1*U(1) + w2*U(2);
dUT_dn1= U(1)*dw1_dn1  + U(2)* dw2_dn1;
dUT_ds1= U(1)*dw1_ds1  + U(2)* dw2_ds1;
dUT_ds2= U(1)*dw1_ds2  + U(2)* dw2_ds2;
dUT_dT1= U(1)*dw1_dT1  + U(2)* dw2_dT1;
dUT_dT2= U(1)*dw1_dT2  + U(2)* dw2_dT2;
dUT_dD1= U(1)*dw1_dD1  + U(2)* dw2_dD1;
dUT_dD2= U(1)*dw1_dD2  + U(2)* dw2_dD2;
dUT_dU1= U(1)*dw1_dU1  + U(2)* dw2_dU1 + w1;
dUT_dU2= U(1)*dw1_dU2  + U(2)* dw2_dU2 + w2;

VT= Vb(1)*w1 + Vb(2)*w2;

%--------------------------- laminar part --------------------------------

% sets up the differences for the laminar part from node 1 to Transition position
hl=sT-s1;
Dl=[D(1); DT];
Tl=[T(1); TT];
Ul=[U(1); UT];
Vl=[Vb(1);VT];
Hl=Dl./Tl;
Rel= Tl.*Ul/nu;

[ f1,f2,der] = JacobiLam( Tl,Ul,Vl,Hl,Rel,hl,s1,nu);


fl=[f1;f2];


% det derivates in respect to 1 and 2 variables

df1_dT1L=der.df1_dT(1) + der.df1_dT(2)*dTT_dT1 + der.df1_dD(2)*dDT_dT1 + der.df1_dU(2)*dUT_dT1 + der.df1_ds(2)*dsT_dT1 ;
df1_dD1L=der.df1_dD(1) + der.df1_dT(2)*dTT_dD1 + der.df1_dD(2)*dDT_dD1 + der.df1_dU(2)*dUT_dD1 + der.df1_ds(2)*dsT_dD1 ;
df1_dU1L=der.df1_dU(1) + der.df1_dT(2)*dTT_dU1 + der.df1_dD(2)*dDT_dU1 + der.df1_dU(2)*dUT_dU1 + der.df1_ds(2)*dsT_dU1 ;
df1_ds1L=der.df1_ds(1) + der.df1_dT(2)*dTT_ds1 + der.df1_dD(2)*dDT_ds1 + der.df1_dU(2)*dUT_ds1 + der.df1_ds(2)*dsT_ds1 ;
df2_dT1L=der.df2_dT(1) + der.df2_dT(2)*dTT_dT1 + der.df2_dD(2)*dDT_dT1 + der.df2_dU(2)*dUT_dT1 + der.df2_ds(2)*dsT_dT1 ;
df2_dD1L=der.df2_dD(1) + der.df2_dT(2)*dTT_dD1 + der.df2_dD(2)*dDT_dD1 + der.df2_dU(2)*dUT_dD1 + der.df2_ds(2)*dsT_dD1 ;
df2_dU1L=der.df2_dU(1) + der.df2_dT(2)*dTT_dU1 + der.df2_dD(2)*dDT_dU1 + der.df2_dU(2)*dUT_dU1 + der.df2_ds(2)*dsT_dU1 ;
df2_ds1L=der.df2_ds(1) + der.df2_dT(2)*dTT_ds1 + der.df2_dD(2)*dDT_ds1 + der.df2_dU(2)*dUT_ds1 + der.df2_ds(2)*dsT_ds1 ;

df1_dT2L=der.df1_dT(2)*dTT_dT2 + der.df1_dD(2)*dDT_dT2 + der.df1_dU(2)*dUT_dT2 + der.df1_ds(2)*dsT_dT2 ;
df1_dD2L=der.df1_dT(2)*dTT_dD2 + der.df1_dD(2)*dDT_dD2 + der.df1_dU(2)*dUT_dD2 + der.df1_ds(2)*dsT_dD2 ;
df1_dU2L=der.df1_dT(2)*dTT_dU2 + der.df1_dD(2)*dDT_dU2 + der.df1_dU(2)*dUT_dU2 + der.df1_ds(2)*dsT_dU2 ;
df1_ds2L=der.df1_dT(2)*dTT_ds2 + der.df1_dD(2)*dDT_ds2 + der.df1_dU(2)*dUT_ds2 + der.df1_ds(2)*dsT_ds2 ;
df2_dT2L=der.df2_dT(2)*dTT_dT2 + der.df2_dD(2)*dDT_dT2 + der.df2_dU(2)*dUT_dT2 + der.df2_ds(2)*dsT_dT2 ;
df2_dD2L=der.df2_dT(2)*dTT_dD2 + der.df2_dD(2)*dDT_dD2 + der.df2_dU(2)*dUT_dD2 + der.df2_ds(2)*dsT_dD2 ;
df2_dU2L=der.df2_dT(2)*dTT_dU2 + der.df2_dD(2)*dDT_dU2 + der.df2_dU(2)*dUT_dU2 + der.df2_ds(2)*dsT_dU2 ;
df2_ds2L=der.df2_dT(2)*dTT_ds2 + der.df2_dD(2)*dDT_ds2 + der.df2_dU(2)*dUT_ds2 + der.df2_ds(2)*dsT_ds2 ;

df1_dn1L=der.df1_dT(2)*dTT_dn1 + der.df1_dD(2)*dDT_dn1 + der.df1_dU(2)*dUT_dn1 + der.df1_ds(2)*dsT_dn1 ;
df2_dn1L=der.df2_dT(2)*dTT_dn1 + der.df2_dD(2)*dDT_dn1 + der.df2_dU(2)*dUT_dn1 + der.df2_ds(2)*dsT_dn1 ;

if IsForced % derivates in respect to forced transition arclength
    dw1_dsF= -1/h ;
    dw2_dsF=  1/h; 
    dTT_dsF=T(1)*dw1_dsF  + T(2)* dw2_dsF;
    dDT_dsF=D(1)*dw1_dsF  + D(2)* dw2_dsF;
    dUT_dsF=U(1)*dw1_dsF  + U(2)* dw2_dsF;

    df1L_dsF= der.df1_dT(2)*dTT_dsF + der.df1_dD(2)*dDT_dsF + der.df1_dU(2)*dUT_dsF + der.df1_ds(2) ;
    df2L_dsF= der.df2_dT(2)*dTT_dsF + der.df2_dD(2)*dDT_dsF + der.df2_dU(2)*dUT_dsF + der.df2_ds(2) ;
end


%--------------------------- turbulent part --------------------------------

% sets up the differences for the turbulent part from Transition position to 2 node
ht=s2-sT;
Dt=[DT ; D(2)];
Tt=[TT ; T(2)];
Ut=[UT ; U(2)];
Vt=[VT ; Vb(2)];
Ht=Dt./Tt;
Ret= Tt.*Ut/nu;

% get Ctau of node T
[CT,dCT_dTT,dCT_dDT,dCT_dUT]= InitialCtau( DT,TT,UT,nu, true );

dCT_dn1= dCT_dTT*dTT_dn1 + dCT_dDT*dDT_dn1 + dCT_dUT*dUT_dn1;
dCT_dT1= dCT_dTT*dTT_dT1 + dCT_dDT*dDT_dT1 + dCT_dUT*dUT_dT1;
dCT_dD1= dCT_dTT*dTT_dD1 + dCT_dDT*dDT_dD1 + dCT_dUT*dUT_dD1;
dCT_dU1= dCT_dTT*dTT_dU1 + dCT_dDT*dDT_dU1 + dCT_dUT*dUT_dU1;
dCT_ds1= dCT_dTT*dTT_ds1 + dCT_dDT*dDT_ds1 + dCT_dUT*dUT_ds1;

dCT_dT2= dCT_dTT*dTT_dT2 + dCT_dDT*dDT_dT2 + dCT_dUT*dUT_dT2;
dCT_dD2= dCT_dTT*dTT_dD2 + dCT_dDT*dDT_dD2 + dCT_dUT*dUT_dD2;
dCT_dU2= dCT_dTT*dTT_dU2 + dCT_dDT*dDT_dU2 + dCT_dUT*dUT_dU2;
dCT_ds2= dCT_dTT*dTT_ds2 + dCT_dDT*dDT_ds2 + dCT_dUT*dUT_ds2;

C=[CT; C2];

[ f1,f2,f3,der ] = JacobiTurb(Dt,Tt,C,Ut,Vt,Ht,Ret,ht,sT,false,false,nu);

ft=[f1;f2;f3];

df1_dT1t=der.df1_dT(1).*dTT_dT1 + der.df1_dD(1).*dDT_dT1 + der.df1_dU(1).*dUT_dT1 + der.df1_ds(1).*dsT_dT1 + der.df1_dC(1).*dCT_dT1;
df1_dD1t=der.df1_dT(1).*dTT_dD1 + der.df1_dD(1).*dDT_dD1 + der.df1_dU(1).*dUT_dD1 + der.df1_ds(1).*dsT_dD1 + der.df1_dC(1).*dCT_dD1;
df1_dU1t=der.df1_dT(1).*dTT_dU1 + der.df1_dD(1).*dDT_dU1 + der.df1_dU(1).*dUT_dU1 + der.df1_ds(1).*dsT_dU1 + der.df1_dC(1).*dCT_dU1;
df1_ds1t=der.df1_dT(1).*dTT_ds1 + der.df1_dD(1).*dDT_ds1 + der.df1_dU(1).*dUT_ds1 + der.df1_ds(1).*dsT_ds1 + der.df1_dC(1).*dCT_ds1;

df2_dT1t=der.df2_dT(1).*dTT_dT1 + der.df2_dD(1).*dDT_dT1 + der.df2_dU(1).*dUT_dT1 + der.df2_ds(1).*dsT_dT1 + der.df2_dC(1).*dCT_dT1;
df2_dD1t=der.df2_dT(1).*dTT_dD1 + der.df2_dD(1).*dDT_dD1 + der.df2_dU(1).*dUT_dD1 + der.df2_ds(1).*dsT_dD1 + der.df2_dC(1).*dCT_dD1;
df2_dU1t=der.df2_dT(1).*dTT_dU1 + der.df2_dD(1).*dDT_dU1 + der.df2_dU(1).*dUT_dU1 + der.df2_ds(1).*dsT_dU1 + der.df2_dC(1).*dCT_dU1;
df2_ds1t=der.df2_dT(1).*dTT_ds1 + der.df2_dD(1).*dDT_ds1 + der.df2_dU(1).*dUT_ds1 + der.df2_ds(1).*dsT_ds1 + der.df2_dC(1).*dCT_ds1;

df3_dT1=der.df3_dT(1).*dTT_dT1 + der.df3_dD(1).*dDT_dT1 + der.df3_dU(1).*dUT_dT1 + der.df3_ds(1).*dsT_dT1 + der.df3_dC(1).*dCT_dT1;
df3_dD1=der.df3_dT(1).*dTT_dD1 + der.df3_dD(1).*dDT_dD1 + der.df3_dU(1).*dUT_dD1 + der.df3_ds(1).*dsT_dD1 + der.df3_dC(1).*dCT_dD1;
df3_dU1=der.df3_dT(1).*dTT_dU1 + der.df3_dD(1).*dDT_dU1 + der.df3_dU(1).*dUT_dU1 + der.df3_ds(1).*dsT_dU1 + der.df3_dC(1).*dCT_dU1;
df3_ds1=der.df3_dT(1).*dTT_ds1 + der.df3_dD(1).*dDT_ds1 + der.df3_dU(1).*dUT_ds1 + der.df3_ds(1).*dsT_ds1 + der.df3_dC(1).*dCT_ds1;

df1_dT2t=der.df1_dT(2) + der.df1_dT(1).*dTT_dT2 + der.df1_dD(1).*dDT_dT2 + der.df1_dU(1).*dUT_dT2 + der.df1_ds(1).*dsT_dT2 + der.df1_dC(1).*dCT_dT2;
df1_dD2t=der.df1_dD(2) + der.df1_dT(1).*dTT_dD2 + der.df1_dD(1).*dDT_dD2 + der.df1_dU(1).*dUT_dD2 + der.df1_ds(1).*dsT_dD2 + der.df1_dC(1).*dCT_dD2;
df1_dU2t=der.df1_dU(2) + der.df1_dT(1).*dTT_dU2 + der.df1_dD(1).*dDT_dU2 + der.df1_dU(1).*dUT_dU2 + der.df1_ds(1).*dsT_dU2 + der.df1_dC(1).*dCT_dU2;
df1_ds2t=der.df1_ds(2) + der.df1_dT(1).*dTT_ds2 + der.df1_dD(1).*dDT_ds2 + der.df1_dU(1).*dUT_ds2 + der.df1_ds(1).*dsT_ds2 + der.df1_dC(1).*dCT_ds2;

df2_dT2t=der.df2_dT(2) + der.df2_dT(1).*dTT_dT2 + der.df2_dD(1).*dDT_dT2 + der.df2_dU(1).*dUT_dT2 + der.df2_ds(1).*dsT_dT2 + der.df2_dC(1).*dCT_dT2;
df2_dD2t=der.df2_dD(2) + der.df2_dT(1).*dTT_dD2 + der.df2_dD(1).*dDT_dD2 + der.df2_dU(1).*dUT_dD2 + der.df2_ds(1).*dsT_dD2 + der.df2_dC(1).*dCT_dD2;
df2_dU2t=der.df2_dU(2) + der.df2_dT(1).*dTT_dU2 + der.df2_dD(1).*dDT_dU2 + der.df2_dU(1).*dUT_dU2 + der.df2_ds(1).*dsT_dU2 + der.df2_dC(1).*dCT_dU2;
df2_ds2t=der.df2_ds(2) + der.df2_dT(1).*dTT_ds2 + der.df2_dD(1).*dDT_ds2 + der.df2_dU(1).*dUT_ds2 + der.df2_ds(1).*dsT_ds2 + der.df2_dC(1).*dCT_ds2;

df3_dT2=der.df3_dT(2) + der.df3_dT(1).*dTT_dT2 + der.df3_dD(1).*dDT_dT2 + der.df3_dU(1).*dUT_dT2 + der.df3_ds(1).*dsT_dT2 + der.df3_dC(1).*dCT_dT2;
df3_dD2=der.df3_dD(2) + der.df3_dT(1).*dTT_dD2 + der.df3_dD(1).*dDT_dD2 + der.df3_dU(1).*dUT_dD2 + der.df3_ds(1).*dsT_dD2 + der.df3_dC(1).*dCT_dD2;
df3_dU2=der.df3_dU(2) + der.df3_dT(1).*dTT_dU2 + der.df3_dD(1).*dDT_dU2 + der.df3_dU(1).*dUT_dU2 + der.df3_ds(1).*dsT_dU2 + der.df3_dC(1).*dCT_dU2;
df3_ds2=der.df3_ds(2) + der.df3_dT(1).*dTT_ds2 + der.df3_dD(1).*dDT_ds2 + der.df3_dU(1).*dUT_ds2 + der.df3_ds(1).*dsT_ds2 + der.df3_dC(1).*dCT_ds2;

df1_dn1t=der.df1_dT(1).*dTT_dn1 + der.df1_dD(1).*dDT_dn1 + der.df1_dU(1).*dUT_dn1 + der.df1_ds(1).*dsT_dn1 + 0;
df2_dn1t=der.df2_dT(1).*dTT_dn1 + der.df2_dD(1).*dDT_dn1 + der.df2_dU(1).*dUT_dn1 + der.df2_ds(1).*dsT_dn1 + der.df2_dC(1).*dCT_dn1;
df3_dn1t=der.df3_dT(1).*dTT_dn1 + der.df3_dD(1).*dDT_dn1 + der.df3_dU(1).*dUT_dn1 + der.df3_ds(1).*dsT_dn1 + der.df3_dC(1).*dCT_dn1;

df2_dC2=der.df2_dC(2);
df3_dC2=der.df3_dC(2);

if IsForced % derivates in respect to forced transition arclength
    dCT_dsF= dCT_dTT*dTT_dsF + dCT_dDT*dDT_dsF + dCT_dUT*dUT_dsF;
    df1t_dsF= der.df1_dT(1)*dTT_dsF + der.df1_dD(1)*dDT_dsF + der.df1_dU(1)*dUT_dsF + der.df1_ds(1) + der.df1_dC(1).*dCT_dsF ;
    df2t_dsF= der.df2_dT(1)*dTT_dsF + der.df2_dD(1)*dDT_dsF + der.df2_dU(1)*dUT_dsF + der.df2_ds(1) + der.df2_dC(1).*dCT_dsF ;
    df3t_dsF= der.df3_dT(1)*dTT_dsF + der.df3_dD(1)*dDT_dsF + der.df3_dU(1)*dUT_dsF + der.df3_ds(1) + der.df3_dC(1).*dCT_dsF ;
    erg.df1_dsF= df1L_dsF + df1t_dsF; 
    erg.df2_dsF= df2L_dsF + df2t_dsF; 
    erg.df3_dsF= df3t_dsF;  
end

% adds turbulent and laminar part of panel equation
% -> Ctau Equation f3 only for turbulent part

f= [fl(1) + ft(1) ;fl(2) + ft(2) ; ft(3) ];

erg.df1_dT= [df1_dT1L + df1_dT1t, df1_dT2L + df1_dT2t]; 
erg.df1_dD= [df1_dD1L + df1_dD1t, df1_dD2L + df1_dD2t]; 
erg.df1_dU= [df1_dU1L + df1_dU1t, df1_dU2L + df1_dU2t]; 
erg.df1_ds= [df1_ds1L + df1_ds1t, df1_ds2L + df1_ds2t]; 

erg.df2_dT= [df2_dT1L + df2_dT1t, df2_dT2L + df2_dT2t]; 
erg.df2_dD= [df2_dD1L + df2_dD1t, df2_dD2L + df2_dD2t]; 
erg.df2_dU= [df2_dU1L + df2_dU1t, df2_dU2L + df2_dU2t]; 
erg.df2_ds= [df2_ds1L + df2_ds1t, df2_ds2L + df2_ds2t]; 

erg.df3_dT= [df3_dT1, df3_dT2]; 
erg.df3_dD= [df3_dD1, df3_dD2]; 
erg.df3_dU= [df3_dU1, df3_dU2]; 
erg.df3_ds= [df3_ds1, df3_ds2]; 

erg.df1_dC= [df1_dn1L + df1_dn1t , 0];
erg.df2_dC= [df2_dn1L + df2_dn1t , df2_dC2];
erg.df3_dC= [df3_dn1t            , df3_dC2];



end

