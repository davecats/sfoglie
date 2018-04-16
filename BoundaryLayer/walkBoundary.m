function sol= walkBoundary(prf,wake,sol,flo,eng,section)
%WALKBOUNDARY integrates the Boundary Layer equations with the pre set initial value and velocity vector 
%             using finite differences method for discretisation and solves the resulting nonlinear equations
%             with a Newton Method for each integration step.
%             Looks for seperation with the Hartree condition (H12 > Hsep, beta<-0.199) 
%             If seperation occurs presets plausible H values and adjust velocity
%             section   1: suctionside
%                       2: pressure side
%                       3: wake
%


% read critical amplifikation factor
nkrit=flo.nkrit;
nu=flo.nu;
Lges=[prf.panels.L, wake.L]; 


% set start + end node index and step direction for each section
% -----------------------------------------------------------------------
if section== 1 % suction side
    lam=true; tr=false; IsWake=false;
    Start=prf.Nle-1; step=-1; Ende=2; shift=1;
    s=prf.sU; 
    sBL=s;% arc length in boundary layer direction
    gap=zeros(prf.Nle-1,1);
elseif section == 2  % pressure side
    lam=true; tr=false;IsWake=false;
    Start=prf.Nle; step=1; Ende=prf.N-1; shift=0;
    s=prf.s; 
    sBL=[zeros(1,prf.Nle-1),prf.sL];
    gap=zeros(prf.N,1);
else  % wake
   lam=false;  IsWake=true;
   Start=prf.N +1; step=1; Ende=prf.N+wake.N-1; shift=0;
   s=wake.s; 
   sBL= [ zeros(1,prf.N),wake.s(1:end-1)+prf.sL(end)]; 
   gap=[zeros(prf.N,1);wake.gap];
   C1= (sol.c(1)*sol.T(1) + sol.c(prf.N)*sol.T(prf.N) ) / (sol.T(1)+sol.T(prf.N)); % initial wake Ctau
   sol.c(prf.N+1)=C1;
   C2=C1;
end


% loop to go through all integration intervalls
% "1" value: known value at intervalls startpoint 
% "2" value: variables to be solved for
for i=Start:step:Ende
    k=0; res=1; % iteration counter and residuum
    % set velocity for 1 and 2 node
    Vb=sol.Vb(i:step:i+step);
    U = sol.U(i:step:i+step);    
    % set values for node 1
    T1=sol.T(i);
    D1=sol.D(i);
    HK1=sol.HK(i);
    % starting values for Newton iteration 2=1
    T2=T1; 
    D2=D1;
    HK2=HK1;
    
    HKset=false;
    
%-----------------------------------------------------LAM----------------------------------------------------------            
    if lam % laminar part of boundary layer
        % Newton method iteration
        while res>sol.resmax && k<sol.itmax
            % updates values
            D=[D1;D2];T=[T1;T2];
            H=[HK1;HK2];
            Ret= T.*U/nu; 
            if HKset==false % no seperation occured
                [ f1,f2,df_dT,df_dD ] = SingleJacobiLam(T,U, Vb,H,Ret,Lges(i-shift),sBL(i),false,nu); 
                % Jacobimatrix for unknown values at 2
                J=[df_dT, df_dD];         

                % correction of old solution x(k+1)=x(k) + dx(k)
                % solution of the Equation J dx(k)=-f(k) with cramers law      
                dT2= (f2*J(1,2) - f1*J(2,2))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );
                dD2= (f1*J(2,1) - f2*J(1,1))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );

                res=max(abs( [dT2/T2,dD2/D2] ));
                % under relaxation for big changes
                if res>0.3; Rel=0.3/res; else Rel=1; end            
                % update values
                T2=T2 + Rel*dT2;
                D2=D2 + Rel*dD2;
            end
            % substrac Trailing edge panel gap for shape parameter
            HK2=( D2-gap(i+step) )/T2;
            H(2)=HK2; 
            
            % check for seperation -> occurs when shape parameter reaches a certain value
            % -> prescribe the growth of H and adapts the tangential boundary edge velocity
            if (HK2 > sol.HmaxLam && HKset==false) ||i==Ende; % force seperation for TE points
                HKset=true;
                Htmp=HK1 + 0.03*Lges(i-shift)/T(1) ; % limit increase of H
                Htmp=max(Htmp,sol.HmaxLam);
                %disp(['seperation at node: ' num2str(i+step)] );
            end 
                

            % seperation 
            %-------------------------------------------------------------------------------
            if HKset
                % add condition H2 = Htmp to Newton System and solves also for new velocity U2
                [ f1,f2,df_dT,df_dD,df_dU] = SingleJacobiLam(T,U, Vb,H,Ret,Lges(i-shift),sBL(i),true,nu); 
                J=[df_dT, df_dD, df_dU; -D(2)/T(2)^2, 1/T(2),0];
                % right hand side of EQ system
                rhs=[-f1;-f2; Htmp - D(2)/T(2)];
                % solves the EQ system
                dz=J\rhs;
                
                dT2=dz(1);dD2=dz(2);dU2=dz(3);
                
                res=max(abs( [dT2/T(2),dD2/D(2),dU2/U(2)] ));
                % update values
                if res>0.3; Rel=0.3/res; else Rel=1; end
                T2=T(2) + Rel*dT2;
                D2=D(2) + Rel*dD2;
                U(2)=U(2)+Rel*dU2;  
            end
            %----------------------------------------------------------------------------------            
            
            % correct for to small H values: adjust D to have H>1.02
            dh= max(0,1.02-( D2-gap(i+step) )/T2);
            D2=D2 +dh*T2;

            k=k+1; %iteration count       
        end % <- End of iteration for current intervall

%------------------------------------------------------TURB---------------------------------------------------------    
    else % turbulent part of boundary layer 
        while res>sol.resmax && k<sol.itmax
            %update values
            D=[D1;D2];T=[T1;T2];
            H=[HK1;HK2];
            Ret= T.*U/nu;
            Ctau=[C1;C2];
            
            if  HKset==false  % no seperation occured
                [ f1,f2,f3,df_dT,df_dD,df_Ct ] = SingleJacobiTurb(D,T,Ctau,U,Vb,H,Ret,Lges(i-shift),sBL(i),false,IsWake,gap(i:i+step),nu); 
                % Jacobimatrix for unknown values at 2
                J=[df_dT, df_dD,df_Ct ];
                rhs=-[f1;f2;f3];
                dz=J\rhs;
                dT2=dz(1);dD2=dz(2); dC2=dz(3);

                res=max(abs( [dT2/T2,dD2/D2,dC2/C2]));
                % under relaxation for big changes
                if res>0.3; Rel=0.3/res; else Rel=1; end
                % update values
                T2=T2 + Rel*dT2;
                D2=D2 + Rel*dD2;
                C2=C2 + Rel*dC2;
            end
            HK2=( D2-gap(i+step) )/T2;
            H(2)=HK2;
            
            
            % check for seperation -> occurs when when shape parameter reaches a certain value
            % -> prescribe th growth of H and adapts the tangential boundary edge velocity
            if HK2 > sol.HmaxTurb && HKset==false; 
                HKset=true;
                if IsWake
                    tmp=0.03*Lges(i-shift)/T(1);
                    Htmp=HK1;
                    Htmp= Htmp - (Htmp + tmp*(Htmp-1)^3 - HK1)/(1+3*tmp*(Htmp-1)^2);
                    Htmp= Htmp - (Htmp + tmp*(Htmp-1)^3 - HK1)/(1+3*tmp*(Htmp-1)^2);
                    Htmp= Htmp - (Htmp + tmp*(Htmp-1)^3 - HK1)/(1+3*tmp*(Htmp-1)^2);
                    Htmp=max(Htmp,1.01);
                else
                    Htmp= HK1-0.15*Lges(i-shift)/T(1); % slow decrease in shape parameter
                    Htmp=max(Htmp,sol.HmaxTurb);
                end
                %disp(['seperation at node: ' num2str(i+step)] );
            end 
                
            
            % seperation 
            %-------------------------------------------------------------------------------
            if HKset
                % add condition H2 = Htmp to Newton System and solves for new velocity U2
                [ f1,f2,f3,df_dT,df_dD,df_Ct,df_dU] = SingleJacobiTurb(D,T,Ctau,U,Vb,H,Ret,Lges(i-shift),sBL(i),true,IsWake,gap(i:i+step),nu);
                % Jacobimatrix for unknown values at 2
                J=[df_dT, df_dD, df_Ct,df_dU; -D(2)/T(2)^2, 1/T(2),0,0];
                rhs=[-f1;-f2;-f3; Htmp- D(2)/T(2)];
                dz=J\rhs;
                
                dT2=dz(1);dD2=dz(2);dC2=dz(3);dU2=dz(4);
                
                res=max(abs( [dT2/T(2),dD2/D(2),dC2/Ctau(2),dU2/U(2)] ));
                if res>0.3; Rel=0.3/res; else Rel=1; end
                T2=T(2) + Rel*dT2;
                D2=D(2) + Rel*dD2;
                C2=Ctau(2) + Rel*dC2;   
                U(2)=U(2) + Rel*dU2;
            end
            %----------------------------------------------------------------------------------
            
            % limit Ctau to realistic values
            C2=min(C2,0.3);
            C2=max(C2,0.0000001);
            
            if IsWake; 
                Hmin=1.00005;
            else
                Hmin=1.02;
            end
            % adjust D to hava H > Hmin if necessery
            dh= max(0,Hmin-( D2-gap(i+step) )/T2);
            D2=D2 + dh*T2;
            
            k=k+1; % iteration counter
        end % <- End of iteration loop for current intervall
        % update Ctau solution
        sol.c(i+step)=C2;C1=C2;
        
    end % End of lam/turb condition
    
%------------------------------------------------- Update panel point 2 solution -------------------------------------       
    % For to big residuals -> extrapolate solution
    if res>0.1 
        disp(['not konverged->solution extrapolated, node: ' num2str(i+step) ', residuum: ' num2str(res)] );
        if IsWake
            T2=T1; 
            tmp=(s(i-prf.N+step)-s(i-prf.N))/(10*D1);
            D2= (D1+T1*tmp)/(1+tmp);
            %D2=D1;
            C1=sol.c(i);C2=C1;
            sol.c(i+step)=C1;
        else
            %T2=T1*sqrt(s(i+step)/s(i));
            Tt=sol.T(i) + Lges(i-shift)*(sol.T(i)-sol.T(i-step))/Lges(i-shift-step);
            T2=Tt;%0.5*(T2+Tt);
            %D2=D1*sqrt(s(i+step)/s(i)); 
            Dt=sol.D(i) + Lges(i-shift)*(sol.D(i)-sol.D(i-step))/Lges(i-shift-step);
            D2=Dt;%0.5*(D2+Dt);
        end
        if  i~=Ende; % central approximation for U
            U(2)= (sol.U(i+2*step)*Lges(i-shift)+sol.U(i)*Lges(i-shift+step))/(Lges(i-shift+step)+Lges(i-shift));
        else % backward approximation for U
            U(2)= sol.U(i);
        end
    end
    %---------------------------------------------
    


    if lam
            % evaluate the amplification equation and check for transition
        %---------------------------------------------------------------------
        n2 = AmplSol(flo,sol.c(i), T, U, H,Ret, Lges(i-shift) );
        sol.c(i+step)=n2;
        
        % tripping arc length in current intervall
        if sol.Tripping(section) && prf.nodes.X(i+step)>sol.xT(section); 
            tr=true;
            w1= (prf.nodes.X(i+step)-sol.xT(section))/ (prf.nodes.X(i+step)-prf.nodes.X(i));
            w2=1-w1;
            sol.sT(section)= w1*sBL(i) + w2*sBL(i+step);
        end

        % transition occures in current intervall
        if sol.c(i+step)> nkrit || tr; 
            lam=false;     
            
            if sol.c(i+step)>nkrit; sol.Tripping(section)=false; end
            
            % solve the Equations for Transition panel
            ind=i:step:i+step; % node indizes
            [ Llam,T2,D2,U(2),C1,~ ] = Transition(section, sol, flo, eng, sol.c(ind), sol.T(ind), sol.D(ind),sol.U(ind),sol.Vb(ind),sBL(i),Lges(i-shift));
            
            %----------
            % write transition node values
            sol.iTran(section)=i+step;
            sol.tran.s(section)=prf.s(i)-Llam;
            sol.tran.n2(section)=sol.c(i+step);
            sol.tran.Llam(section)=Llam;
            sol.tran.Lturb(section)=Lges(i-shift)-Llam;
            
            % store Ctau instead of n 
            C2=C1;
            sol.c(i+step)=C1; 
        end

       %if no transition occured so far -> set transition to last airfoil node
       if i==Ende
            sol.Tripping(section)=true;
            lam=false;
            sol.sT(section)=sBL(i+step);
            sol.xT(section)=prf.nodes.X(i+step);
            
            % set initial value for sqrt(Ctau) at TE
            C1= InitialCtau( D(2),T(2),U(2),nu );
            
            %----------            
            % write transition node values
            sol.iTran(section)=i+step;
            sol.tran.s(section)=prf.s(i+step);
            sol.tran.n2(section)=sol.c(i+step);
            sol.tran.Llam(section)=Lges(i-shift);
            sol.tran.Lturb(section)=0;
            
            % store Ctau instead of n 
            C2=C1;
            sol.c(i+step)=C1; 
       end

    end
    %----------------------------------------------------------------------------------------------------------------------------------
    sol.T(i+step)=T2;
    sol.D(i+step)=D2;
    sol.U(i+step)=U(2);
    sol.HK(i+step)=( D2-gap(i+step) )/T2;
    sol.Ret(i+step)=T2*U(2)/nu;
end  % <- End of the loop over all boundary subintervalls


end






