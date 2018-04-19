function [ f, erg,sT ] = TransitionEQ2( flo,eng, n, T, D,U,Vb,s1,h, C2,sF , IsForced)
%TRANSITIONEQ   gives the equationsystem for the transition panel
%                   different method for determination of transition location
%               -> uses sum of laminar and turbulent part for momentum and shapeparameter equation
%               -> uses Ctau Equation only for turbulent part

nkrit=flo.nkrit;
nu=flo.nu;

if nargin==10; IsForced=false;  end

s2=s1 + h;



%===============================================================================================
%   
%===============================================================================================

if IsForced
    % set weightening factors in case of forced transition
    sT= sF;
else
    
    % secant law approximation for Transitionpanel
     k=0;res=1;
    w1= (n(2)-nkrit)/(n(2)-n(1));
    w2= 1 - w1;
    TT=w1*T(1)+w2*T(2);
    DT=w1*D(1)+w2*D(2);
    UT=w1*U(1)+w2*U(2);
    sT=w1*s1 + w2*s2;
    nT=nkrit;
    hL=sT-s1;
    
    % find exakt point of Transition where n2 = nkrit;
    while res>5e-5 && k<30
        
        V = [Vb(1); w1*Vb(1)+w2*Vb(2)];
        TL=[T(1);TT];
        DL=[D(1);DT];
        UL=[U(1);UT];
        nL=[n(1);nT];
        HL=DL./TL;
        ReL=UL.*TL/nu;
        
    
        [ f1,f2,der] = JacobiLam( TL,UL,V,HL,ReL,hL,s1,nu);
    
        [ f3, df3_dn,df3_dT,df3_dD,~, df3_ds ]=AmplificationEquation(flo,nL,TL,UL,HL,ReL,hL);
        
        J=[der.df1_dT(2),der.df1_dD(2),0        ,der.df1_ds(2);...  % momentum EQ
           der.df2_dT(2),der.df2_dD(2),0        ,der.df2_ds(2);...  % kin Energy EQ
           df3_dT(2)    ,df3_dD(2)    ,df3_dn(2),df3_ds(2) ;...     % Ampl EQ
           0            ,0            ,1     ,0];% force n2=nkrit
        rhs=-[f1;f2;f3;nT-nkrit];
        
        dz=J\rhs;
        
        res=max(abs( [dz(1)/TT,dz(2)/DT,dz(3)/nT,dz(4)/hL] ));
        % under relaxation for big changes
        if res>0.3; rel=0.3/res; else rel=1; end  
        TT = TL(2) + rel*dz(1);
        DT = DL(2) + rel*dz(2);
        nT = nL(2) + rel*dz(3);
        hL = hL    + rel*dz(4);
        
        w2= hL /h;
        w1=1-w2;
        
        k=k+1;
    end
    
    if hL<h && hL>0 % only overwrite approximation if solution makes sense
       sT=s1 + hL;
    else % switch to sekant Method
%------------------------------------------------------------------

        hL=sT-s1;hLmin=hL;
        n1=n(1);n2=n(2);
        h1=0;h2=h;
        k=0;res=1;resmin=2;
        while res>5e-5 && k<20
            w2= hL/h;
            w1= 1 - w2;
            TT=w1*T(1)+w2*T(2);
            DT=w1*D(1)+w2*D(2);
            UT=w1*U(1)+w2*U(2);
            nT=AmplSol(flo,n(1),[T(1);TT], [U(1); UT],[D(1)/T(1); DT/TT],[U(1)*T(1);UT*TT]/nu,hL);
            res=abs(nT-nkrit); 
            if res>resmin; hL=hLmin; break; else resmin=res;hLmin=hL; end
            if nT > nkrit + 5e-5
                h2=hL;n2=nT;
            elseif nT < nkrit - 5e-5
                h1=hL;n1=nT;
            else
               break; 
            end
            %korrektion
            a= (n2-n1)/(h2-h1);
            hL=(nkrit-n2)/a+h2;
            hLmin=hL;
            k=k+1;
        end
        if hL<h && hL>0;
            sT=s1 + hLmin; 
        end;
 %------------------------------------------------------------------       
    end

end


%------------------------------------------------------------------------------------------------------
%                       get the boundary layer equations for Transition panel



% Transition point value derivates
w2= (sT-s1) /h;
w1=1-w2;

if w2<1e-8; w1=1; w2=0; sT=s1; end
if w1<1e-8; w1=0; w2=1; sT=s2; end

%dw2_dsT=1/h;

dw2_ds1=  (w2-1) /h;
dw2_ds2=  -w2    /h;

dw1_ds1= -dw2_ds1;
dw1_ds2= -dw2_ds2;


TT=w1*T(1) + w2*T(2);
dTT_dT1=  w1;
dTT_dT2=  w2;
dTT_ds1=  T(1)*dw1_ds1 +  T(2)*dw2_ds1;
dTT_ds2=  T(1)*dw1_ds2 +  T(2)*dw2_ds2;

DT=w1*D(1) + w2*D(2);
dDT_dD1=  w1;
dDT_dD2=  w2;
dDT_ds1=  D(1)*dw1_ds1 +  D(2)*dw2_ds1;
dDT_ds2=  D(1)*dw1_ds2 +  D(2)*dw2_ds2;


UT=w1*U(1) + w2*U(2);
dUT_dU1= w1;
dUT_dU2= w2;
dUT_ds1=  U(1)*dw1_ds1 +  U(2)*dw2_ds1;
dUT_ds2=  U(1)*dw1_ds2 +  U(2)*dw2_ds2;


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

df1_dT1L=der.df1_dT(1) ;
df1_dD1L=der.df1_dD(1) ;
df1_dU1L=der.df1_dU(1) ;
df1_ds1L=der.df1_ds(1) + der.df1_dT(2)*dTT_ds1 + der.df1_dD(2)*dDT_ds1 + der.df1_dU(2)*dUT_ds1  ;
df2_dT1L=der.df2_dT(1) ;
df2_dD1L=der.df2_dD(1) ;
df2_dU1L=der.df2_dU(1) ;
df2_ds1L=der.df2_ds(1) + der.df2_dT(2)*dTT_ds1 + der.df2_dD(2)*dDT_ds1 + der.df2_dU(2)*dUT_ds1  ;

df1_dT2L=der.df1_dT(2)*dTT_dT2 ;
df1_dD2L=der.df1_dD(2)*dDT_dD2 ;
df1_dU2L=der.df1_dU(2)*dUT_dU2 ;
df1_ds2L=der.df1_dT(2)*dTT_ds2 + der.df1_dD(2)*dDT_ds2 + der.df1_dU(2)*dUT_ds2 ;
df2_dT2L=der.df2_dT(2)*dTT_dT2 ;
df2_dD2L=der.df2_dD(2)*dDT_dD2 ;
df2_dU2L=der.df2_dU(2)*dUT_dU2 ;
df2_ds2L=der.df2_ds(2) + der.df2_dT(2)*dTT_ds2 + der.df2_dD(2)*dDT_ds2 + der.df2_dU(2)*dUT_ds2 ;

df1_dn1L=0 ;
df2_dn1L=0 ;


dw1_dsF= -1/h ;
dw2_dsF=  1/h; 
dTT_dsF=T(1)*dw1_dsF  + T(2)* dw2_dsF;
dDT_dsF=D(1)*dw1_dsF  + D(2)* dw2_dsF;
dUT_dsF=U(1)*dw1_dsF  + U(2)* dw2_dsF;

df1L_dsF= der.df1_dT(2)*dTT_dsF + der.df1_dD(2)*dDT_dsF + der.df1_dU(2)*dUT_dsF + der.df1_ds(2) ;
df2L_dsF= der.df2_dT(2)*dTT_dsF + der.df2_dD(2)*dDT_dsF + der.df2_dU(2)*dUT_dsF + der.df2_ds(2) ;



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


dCT_dT1= dCT_dTT*dTT_dT1;
dCT_dD1= dCT_dDT*dDT_dD1;
dCT_dU1= dCT_dUT*dUT_dU1;
dCT_ds1= dCT_dTT*dTT_ds1 + dCT_dDT*dDT_ds1 + dCT_dUT*dUT_ds1;

dCT_dT2= dCT_dTT*dTT_dT2 ;
dCT_dD2= dCT_dDT*dDT_dD2 ;
dCT_dU2= dCT_dUT*dUT_dU2 ;
dCT_ds2= dCT_dTT*dTT_ds2 + dCT_dDT*dDT_ds2 + dCT_dUT*dUT_ds2;

C=[CT; C2];

[ f1,f2,f3,der ] = JacobiTurb(Dt,Tt,C,Ut,Vt,Ht,Ret,ht,sT,false,false,nu);

ft=[f1;f2;f3];

df1_dT1t= der.df1_dT(1).*dTT_dT1 + der.df1_dC(1).*dCT_dT1;
df1_dD1t= der.df1_dD(1).*dDT_dD1 + der.df1_dC(1).*dCT_dD1;
df1_dU1t= der.df1_dU(1).*dUT_dU1 + der.df1_dC(1).*dCT_dU1;
df1_ds1t= der.df1_dT(1).*dTT_ds1 + der.df1_dD(1).*dDT_ds1 + der.df1_dU(1).*dUT_ds1  + der.df1_dC(1).*dCT_ds1;


df2_dT1t= der.df2_dT(1).*dTT_dT1 + der.df2_dC(1).*dCT_dT1;
df2_dD1t= der.df2_dD(1).*dDT_dD1 + der.df2_dC(1).*dCT_dD1;
df2_dU1t= der.df2_dU(1).*dUT_dU1 + der.df2_dC(1).*dCT_dU1;
df2_ds1t= der.df2_dT(1).*dTT_ds1 + der.df2_dD(1).*dDT_ds1 + der.df2_dU(1).*dUT_ds1  + der.df2_dC(1).*dCT_ds1;

df3_dT1= der.df3_dT(1).*dTT_dT1 + der.df3_dC(1).*dCT_dT1;
df3_dD1= der.df3_dD(1).*dDT_dD1 + der.df3_dC(1).*dCT_dD1;
df3_dU1= der.df3_dU(1).*dUT_dU1 + der.df3_dC(1).*dCT_dU1;
df3_ds1= der.df3_dT(1).*dTT_ds1 + der.df3_dD(1).*dDT_ds1 + der.df3_dU(1).*dUT_ds1  + der.df3_dC(1).*dCT_ds1;

df1_dT2t= der.df1_dT(2) + der.df1_dT(1).*dTT_dT2 + der.df1_dC(1).*dCT_dT2;
df1_dD2t= der.df1_dD(2) + der.df1_dD(1).*dDT_dD2 + der.df1_dC(1).*dCT_dD2;
df1_dU2t= der.df1_dU(2) + der.df1_dU(1).*dUT_dU2 + der.df1_dC(1).*dCT_dU2;
df1_ds2t= der.df1_ds(2) + der.df1_dT(1).*dTT_ds2 + der.df1_dD(1).*dDT_ds2 + der.df1_dU(1).*dUT_ds2 + der.df1_dC(1).*dCT_ds2;

df2_dT2t= der.df2_dT(2) + der.df2_dT(1).*dTT_dT2 + der.df2_dC(1).*dCT_dT2;
df2_dD2t= der.df2_dD(2) + der.df2_dD(1).*dDT_dD2 + der.df2_dC(1).*dCT_dD2;
df2_dU2t= der.df2_dU(2) + der.df2_dU(1).*dUT_dU2 + der.df2_dC(1).*dCT_dU2;
df2_ds2t= der.df2_ds(2) + der.df2_dT(1).*dTT_ds2 + der.df2_dD(1).*dDT_ds2 + der.df2_dU(1).*dUT_ds2 + der.df2_dC(1).*dCT_ds2;

df3_dT2= der.df3_dT(2) + der.df3_dT(1).*dTT_dT2 + der.df3_dC(1).*dCT_dT2;
df3_dD2= der.df3_dD(2) + der.df3_dD(1).*dDT_dD2 + der.df3_dC(1).*dCT_dD2;
df3_dU2= der.df3_dU(2) + der.df3_dU(1).*dUT_dU2 + der.df3_dC(1).*dCT_dU2;
df3_ds2= der.df3_ds(2) + der.df3_dT(1).*dTT_ds2 + der.df3_dD(1).*dDT_ds2 + der.df3_dU(1).*dUT_ds2 + der.df3_dC(1).*dCT_ds2;

df1_dn1t=0;
df2_dn1t=0;
df3_dn1t=0;

df2_dC2=der.df2_dC(2);
df3_dC2=der.df3_dC(2);


dCT_dsF= dCT_dTT*dTT_dsF + dCT_dDT*dDT_dsF + dCT_dUT*dUT_dsF;
df1t_dsF= der.df1_dT(1)*dTT_dsF + der.df1_dD(1)*dDT_dsF + der.df1_dU(1)*dUT_dsF + der.df1_ds(1) + der.df1_dC(1).*dCT_dsF ;
df2t_dsF= der.df2_dT(1)*dTT_dsF + der.df2_dD(1)*dDT_dsF + der.df2_dU(1)*dUT_dsF + der.df2_ds(1) + der.df2_dC(1).*dCT_dsF ;
df3t_dsF= der.df3_dT(1)*dTT_dsF + der.df3_dD(1)*dDT_dsF + der.df3_dU(1)*dUT_dsF + der.df3_ds(1) + der.df3_dC(1).*dCT_dsF ;
erg.df1_dsF= df1L_dsF + df1t_dsF; 
erg.df2_dsF= df2L_dsF + df2t_dsF; 
erg.df3_dsF= df3t_dsF;  


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

