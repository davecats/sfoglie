function [ Llam, T2,D2,U2,C2, der] = RefreshTransition(sec, sol,n, T, D,U,Vb,s1,h,C2 )
%TRANSITION handles the Refresh of the Transition intervall

nu=evalin('base','nu');

C2=min(C2,0.3);
C2=max(C2,0.0000001);

D2=D(2);
T2=T(2);

H2=D(2)/T(2);
Href=H2;
U2=U(2);
Uref=U2;


s2=s1+h;
if sol.Tripping(sec); 
    tr=true; 
else
    tr=false; 
end


% iterate over transition panel, -> Equation is sum of laminar and turbulent part of the transition panel
% f1(h) = f1l(hLam) + f1t(hTurb)
% f2(h) = f2l(hLam) + f2t(hTurb)
% f3(h) = f3t(hTurb) -> Ctau EQ only for turbulent part

k=0; res=1;
 while res>sol.resmax && k<sol.itmax
     
    D(2)=D2; T(2)=T2; U(2)=U2; 
    if tr % Tripping
        [f, der, sT]=TransitionEQ( n, T, D,U,Vb,s1,h, C2, sol.sT(sec),true );
    else
        [f, der, sT]=TransitionEQ( n, T, D,U,Vb,s1,h, C2);
    end
    J= [ der.df1_dT(2), der.df1_dD(2), 0            ,der.df1_dU(2) ; ...
         der.df2_dT(2), der.df2_dD(2), der.df2_dC(2),der.df2_dU(2) ; ...
         der.df3_dT(2), der.df3_dD(2), der.df3_dC(2),der.df3_dU(2) ; ...
         -H2./T2      , 1./T2        , 0            ,0 ];

    rhs=[-f; 1];
    d=J\rhs;
                
    if k<16 % sensitivity to change of H
        senTMP= 1000 * d(4)* Href/Uref;
        if k<6; sen=senTMP; else sen= 0.5*(sen+senTMP); end
    end
    J(4,1)=J(4,1)*Href;
    J(4,2)=J(4,2)*Href; 
    J(4,4)=sen/Uref;
    rhs(4)=-Href^2*(H2/Href -1) - sen*(U2/Uref -1);

    dz=J\rhs;

    res=max(abs( [dz(1)/T2,dz(2)/D2,dz(3)/C2,dz(4)/U2] ));
    if res>0.3; Rel=0.3/res; else Rel=1; end
    T2=T2 + Rel*dz(1);
    D2=D2 + Rel*dz(2);
    C2=C2 + Rel*dz(3);   
    U2=U2 + Rel*dz(4);  

    dh= max(0,1.02-D2/T2);
    D2=D2 +dh*T2;
    C2=min(C2,0.3);
    C2=max(C2,0.0000001);
    
    H2= D2/T2; 
    
    
    k=k+1;      
end % <- End of iteration for current intervall

if res>0.1 
    %disp(['not konverged->solution extrapolated at transition node residuum: ' num2str(res)] );
    T2=T(1)*sqrt(s2/s1);
    D2=D(1)*sqrt(s2/s1); 
    C2=0.05;
    U2=U(2);
end      

Llam= sT-s1;     

 

end

