function sol = PlateINT( sol, x, dx, nu )
%PLATEINT integrates the Boundary-Layer equations for a flat plate and turbulent flow



for i=1:length(sol.T)-1
    HKset=false;
    k=0; res=1; % iteration counter and residuum

    T1=sol.T(i);
    D1=sol.D(i);
    H1=sol.HK(i);
    C1=sol.c(i);
    % starting values for Newton iteration 2=1
    T2=T1; 
    D2=D1;
    H2=H1;
    C2=C1;

    Vb=sol.Vb(i:i+1);
    U=sol.U(i:i+1);
    while res>1e-12 && k<25
        % updates values
        D=[D1;D2];T=[T1;T2];
        H=[H1;H2];C=[C1;C2];
        Ret= T.*U/nu;

        z=[T2;D2;C2];

        if  HKset==false 
            [ f1,f2,f3,df_dT,df_dD,df_Ct ] = SingleJacobiTurb(D,T,C,U,Vb,H,Ret,dx(i),x(i),false,false,[0;0]); 
            % Jacobimatrix for unknown values at 2
            J=[df_dT, df_dD,df_Ct ];
            rhs=-[f1;f2;f3];
            dz=J\rhs;



            res=max(abs(dz./z));
            % under relaxation for big changes
            if res>0.3; Rel=0.3/res; else Rel=1; end

            T2 = z(1) + Rel*dz(1);
            D2 = z(2) + Rel*dz(2);
            C2 = z(3) + Rel*dz(3);
        end
        H2=D2/T2; 

      
            % check for seperation -> occurs when when shape parameter reaches a certain value
            % -> prescribe th growth of H and adapts the tangential boundary edge velocity
            if H2 > sol.HmaxTurb && HKset==false; 
                HKset=true;
                Htmp= H1-0.15*dx(i)/T(1); % slow decrease in shape parameter
                Htmp=max(Htmp,sol.HmaxTurb);
                %disp(['seperation at node: ' num2str(i+step)] );
            end 
            % seperation 
            %-------------------------------------------------------------------------------
            if HKset
                % add condition H2 = Htmp to Newton System and solves for new velocity U2
                [ f1,f2,f3,df_dT,df_dD,df_Ct,df_dU] = SingleJacobiTurb(D,T,C,U,Vb,H,Ret,dx(i),x(i),true,false,[0;0]);
                % Jacobimatrix for unknown values at 2
                J=[df_dT, df_dD, df_Ct,df_dU; -D(2)/T(2)^2, 1/T(2),0,0];
                rhs=[-f1;-f2;-f3; Htmp- D(2)/T(2)];
                dz=J\rhs;
                
                dT2=dz(1);dD2=dz(2);dC2=dz(3);dU2=dz(4);
                
                res=max(abs( [dT2/T(2),dD2/D(2),dC2/C(2),dU2/U(2)] ));
                if res>0.3; Rel=0.3/res; else Rel=1; end
                T2=T(2) + Rel*dT2;
                D2=D(2) + Rel*dD2;
                C2=C(2) + Rel*dC2;   
                U(2)=U(2) + Rel*dU2;
            end
            % limit Ctau to realistic values
            C2=min(C2,0.3);
            C2=max(C2,0.0000001);

            Hmin=1.02;
          
            % adjust D to hava H > Hmin if necessery
            dh= max(0,Hmin-D2/T2);
            D2=D2 + dh*T2;
          k=k+1;
    end
    
        % For to big residuals -> extrapolate solution
    if res>0.1 
        disp(['not konverged->solution extrapolated, node: ' num2str(i+1) ', residuum: ' num2str(res)] );
        %T2=T1*sqrt(s(i+step)/s(i));
        Tt=sol.T(i) + (sol.T(i)-sol.T(i-1));
        T2=Tt;%0.5*(T2+Tt);
        %D2=D1*sqrt(s(i+step)/s(i)); 
        Dt=sol.D(i) + (sol.D(i)-sol.D(i-1));
        D2=Dt;%0.5*(D2+Dt);         
    end
  %----------------------------------------------  
    
    sol.c(i+1)=C2; 
    sol.D(i+1)=D2;  
    sol.T(i+1)=T2; 
    sol.HK(i+1)=D2/T2; 
    sol.Ret(i+1)=U(2)*T2/nu; 
    
end



end

