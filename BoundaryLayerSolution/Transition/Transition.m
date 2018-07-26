function [ Llam, T2,D2,U2,C2, der] = Transition(sec, sol,flo, eng, n, T, D,U,Vb,s1,h)
%TRANSITION Integrates the governing equations over the transition panel.
%The transition point is determined and the panel is split in a laminar and
%a turbulent part.


D2=D(1);
T2=T(1);

s2=s1 + h;

if sol.Tripping(sec); tr=true; else  tr=false; end

% initial value for Ctau at 2 node
C2t= 0.03;

HK1= D(1)/T(1);
HKset=false;

% iterate over transition panel, -> Equation is sum of laminar and turbulent part of the transition panel
% f1(h) = f1l(hLam) + f1t(hTurb)
% f2(h) = f2l(hLam) + f2t(hTurb)
% f3(h) = f3t(hTurb) -> Ctau EQ only for turbulent part

k=0; res=1;
 while res>sol.resmax && k<sol.itmax
     
    D(2)=D2; T(2)=T2; C2=C2t;
    if tr % Tripping
        [f, der, sT]=TransitionEQ( flo, eng, n, T, D,U,Vb,s1,h, C2, sol.sT(sec),true );
    else
        [f, der, sT]=TransitionEQ( flo, eng, n, T, D,U,Vb,s1,h, C2 );
    end
    if HKset==false % no seperation occured
        
        J= [ der.df1_dT(2), der.df1_dD(2), 0 ; ...
             der.df2_dT(2), der.df2_dD(2), der.df2_dC(2) ; ...
             der.df3_dT(2), der.df3_dD(2), der.df3_dC(2) ; ...
            ];
        
        rhs=-f;
        dz=J\rhs;
        res=max(abs( [dz(1)/T2,dz(2)/D2,dz(3)/C2]));
        % under relaxation for big changes
        if res>0.3; Rel=0.3/res; else Rel=1; end
        % update values
        T2=T2 + Rel*dz(1);
        D2=D2 + Rel*dz(2);
        C2t=C2+ Rel*dz(3);
    end
    HK2=D2/T2;

    % check for seperation -> occurs when when shape parameter reaches a certain value
    % -> prescribe th growth of H and adapts the tangential boundary edge velocity
    if HK2 > sol.HmaxLam && HKset==false; 
        HKset=true;
        Htmp=HK1 + 0.03*(sT-s1)/T(1) -0.15*(s2-sT)/T(1); % limit increase of H
        Htmp=max(Htmp,sol.HmaxLam);
    end 

    % seperation 
    %-------------------------------------------------------------------------------
    if HKset
         J= [ der.df1_dT(2), der.df1_dD(2), 0            , der.df1_dU(2); ...
              der.df2_dT(2), der.df2_dD(2), der.df2_dC(2), der.df2_dU(2); ...
              der.df3_dT(2), der.df3_dD(2), der.df3_dC(2), der.df3_dU(2); ...
              -D(2)/T(2)^2 , 1/T(2)       ,0             ,0 ];
            
        
        rhs=[-f;Htmp- D(2)/T(2)];
        dz=J\rhs;
        res=max(abs( [dz(1)/T2,dz(2)/D2,dz(3)/C2, dz(4)/U(2)])); 
        % under relaxation for big changes
        if res>0.3; Rel=0.3/res; else Rel=1; end
        % update values
        T2=T2 + Rel*dz(1);
        D2=D2 + Rel*dz(2);
        C2t=C2 + Rel*dz(3);
        U(2)=U(2) + Rel*dz(4); 
    end
    %----------------------------------------------------------------------------------            
    
    dh= max(0,1.02-D2/T2);
    D2=D2 +dh*T2;

    
    k=k+1;      
end % <- End of iteration for current intervall


Llam= sT-s1;     

U2=U(2);
 

end

