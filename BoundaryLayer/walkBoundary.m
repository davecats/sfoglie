function sol= walkBoundary(prf,wake,sol,section)
%WALKBOUNDARY integrates the Boundary Layer equations with the pre set initial value using finite 
%             differences method for discretisation and solves the resulting nonlinear equations
%             with a Newton Method for each integration step
%             section   1: suctionside
%                       2: pressure side
%                       3: wake
%

nu=evalin('base','nu');
Lges=[prf.panels.L, wake.L]; 
LTE= [zeros(size(prf.panels.L)) , wake.GAP];

% set start + end node index and step direction for each section
% -----------------------------------------------------------------------
if section== 1 % suction side
    lam=true; IsWake=false;
    Start=prf.Nle-1; step=-1; Ende=2; shift=1;
    s=prf.so; 
    
elseif section == 2  % pressure side
    lam=true; IsWake=false;
    Start=prf.Nle; step=1; Ende=prf.N-1; shift=0;
    s=prf.s; 
    
else  % wake
   lam=false;  IsWake=true;
   Start=prf.N +1; step=1; Ende=prf.N+wake.N-1; shift=0;
   s=wake.s; 
   C1= (sol.c(1)*sol.T(1) + sol.c(prf.N)*sol.T(prf.N) ) / (sol.T(1)+sol.T(prf.N)); % initial wake Ctau
   sol.c(prf.N+1)=C1;
   C2=C1;
end


% loop to go through all integration intervalls
for i=Start:step:Ende
    k=0; res=1; % iteration counter and residuum
    % set velocity for 1 and 2 node
    Vb1=sol.Vb(i);Vb2=sol.Vb(i+step);
    Vb=[Vb1;Vb2];
    U1=sol.U(i);U2=sol.U(i+step);
    U=[U1;U2];
    % set values for node 1
    T1=sol.T(i);
    T2=T1; % starting for Newton iteration 2=1
    D1=sol.D(i);
    D2=D1;
    
    HK1=sol.HK(i);
    
    HKset=false;
    
%-----------------------------------------------------LAM----------------------------------------------------------            
    if lam % laminar part of boundary layer
        while res>sol.resmax && k<sol.itmax
            % updates values
            D=[D1;D2];T=[T1;T2];
            
            %if HKset==false % no seperation occured

                [ f1,f2,df_dT,df_dD, params ] = SingleJacobiLam(D,T,U, Vb,Lges(i-shift),[LTE(i); LTE(i+step)],false,IsWake); 
                % Jacobimatrix for unknown values at 2
                J=[df_dT, df_dD];

                % correction of old solution x(k+1)=x(k) + dx(k)
                % solution of the Equation J dx(k)=-f(k) with cramers law      
                dT2= (f2*J(1,2) - f1*J(2,2))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );
                dD2= (f1*J(2,1) - f2*J(1,1))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );

                res=max(abs(dT2/T2),abs(dD2/D2));
                % under relaxation for big changes
                if res>0.3; Rel=0.3/res; else Rel=1; end            
                % update values
                T2=T2 + Rel*dT2;
                D2=D2 + Rel*dD2;

                HK2=D2/T2;
                sol.HK(i+step)=HK2;
                
                % check for seperation -> occurs when when shape parameter reaches a certain value
                % -> prescribe th growth of H and adapts the tangential boundary edge velocity
                if HK2 > sol.HmaxLam && HKset==false; 
                    HKset=true;
                    Htmp=HK1 + 0.03*Lges(i-shift)/T(1) ; % limit increase of H
                    Htmp=max(Htmp,sol.HmaxLam);
                    disp(['HK set at node ' num2str(i+step) ' HK=' num2str(Htmp)  ] );
                end 
                
            %end
            
            % seperation 
            %-------------------------------------------------------------------------------
            if HKset

                % add condition H2 = Htmp to Newton System and solves for new velocity U2
                [ f1,f2,df_dT,df_dD, params ,df_dU] = SingleJacobiLam(D,T,U, Vb,Lges(i-shift),[LTE(i); LTE(i+step)],true,IsWake);   
                
                J=[df_dT, df_dD, df_dU; -D(2)/T(2)^2, 1/T(2),0];
                rhs=[-f1;-f2; Htmp - D(2)/T(2)];
                dz=J\rhs;
                
                dT2=dz(1);dD2=dz(2);dU2=dz(3);
                res=max(abs(dT2/T(2)),abs(dD2/D(2)));
                res=max(res,abs(dU2/U(2)));
                % update values
                if res>0.3; Rel=0.3/res; else Rel=1; end
                T2=T(2)+Rel*dT2;
                D2=D(2)+Rel*dD2;
                U2=U(2)+Rel*dU2;
                U(2)=U2;               
            end
            %----------------------------------------------------------------------------------            
            
            % correct for to small H values
            dh= max(0,1.02-D2/T2);
            D2=D2 +dh*T2;

            k=k+1; %iteration count       
        end % <- End of iteration for current intervall

        
        % evaluate the amplification equation and check for transition
        %---------------------------------------------------------------------
        sol.c(i+step)=IntegrateAmplification(sol.c(i),Lges(i-shift),T1,T2,D1/T1,D2/T2,U1*T1/nu,U2*T2/nu); 
        
        % transition occures in current intervall
        if sol.c(i+step)>9; 
            Cstart=0.03; % set initial value for sqrt(Ctau)
            % split intervall in laminar and turbulent part and calculate new 2 values
            [ Ttran,T2,Dtran,D2,U2,Utran,C1,Llam] =Transition(sol,i,step,Cstart,U,Vb,Lges(i-shift));
            lam=false; 
            
            % Testwerte XFoil
            %D2=2.313544223648736e-002;
            %T2=2.344335248645479e-003;
            %U2=1.28078788668137;
            %C1=0.149801814739205;

            sol.HK(i+step)=D2/T2;
            
            %----------
            % write transition node values
            if section==1;
                sol.TranU=i+step; sol.tran.sU= prf.s(i)-Llam; sol.tran.TU=Ttran;sol.tran.DU=Dtran;sol.tran.UU=Utran;
            else
                sol.TranL=i+step; sol.tran.sL= prf.s(i)+Llam; sol.tran.TL=Ttran;sol.tran.DL=Dtran;sol.tran.UL=Utran;
            end 
            C2=C1;
            sol.c(i+step)=C1; 
        end

        %if no transition occured so far -> set transition to last airfoil node
        if i+step==1 || i+step==prf.N
            % set initial value for sqrt(Ctau)
            fac=1.8*exp(-3.3./(D./T-1));
            [ HS, ~, ~  ]=H32turb( D./T,U.*T/nu);
            Us=0.5*HS.*( -1/3 + T./(0.75*D) );
            [CEQ,~,~,~,~]=CtEQ( D./T,U.*T/nu,HS,Us,false);
            C1=fac(1)* CEQ(1);
            lam=false;
            if section==1; sol.TranU=i+step; else sol.TranL=i+step;  end 
            C2=C1;
            sol.c(i+step)=C1; 
        end 
%------------------------------------------------------TURB---------------------------------------------------------    
    else % turbulent part of boundary layer 
%         if i==9
%            debug=true; 
%         end
        
        while res>sol.resmax && k<sol.itmax
            D=[D1;D2];T=[T1;T2];
            Ctau=[C1;C2];
            
            %if HKset==false % no seperation occured

                [ f1,f2,f3,df_dT,df_dD,df_Ct, params ] = SingleJacobiTurb(D,T,Ctau,U,Vb,Lges(i-shift),[LTE(i); LTE(i+step)],false,IsWake); 
                % Jacobimatrix for unknown values at 2
                J=[df_dT, df_dD,df_Ct ];
                rhs=-[f1;f2;f3];
                dz=J\rhs;
                dT2=dz(1);dD2=dz(2); dC2=dz(3);

                res=max(abs(dT2/T2),abs(dD2/D2));
                res=max(res,abs(dC2/C2));
                % under relaxation for big changes
                if res>0.3; Rel=0.3/res; else Rel=1; end
                % update values
                T2=T2+Rel*dT2;
                D2=D2+Rel*dD2;
                C2=C2+Rel*dC2;

                HK2=D2/T2;
                sol.HK(i+step)=HK2;
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
                    disp(['HK set at node ' num2str(i+step) ' HK=' num2str(Htmp)  ] ); 
                end 
                
            %end
            
            % seperation 
            %-------------------------------------------------------------------------------
            if HKset%HK2 > sol.HmaxTurb%
                % add condition H2 = Htmp to Newton System and solves for new velocity U2
                [ f1,f2,f3,df_dT,df_dD,df_Ct, params,df_dU] = SingleJacobiTurb(D,T,Ctau,U,Vb,Lges(i-shift),[LTE(i); LTE(i+step)],true,IsWake);
                
                J=[df_dT, df_dD, df_Ct,df_dU; -D(2)/T(2)^2, 1/T(2),0,0];
                rhs=[-f1;-f2;-f3; Htmp- D(2)/T(2)];
                dz=J\rhs;
                dT2=dz(1);dD2=dz(2);dC2=dz(3);dU2=dz(4);
                res=max(abs(dT2/T(2)),abs(dD2/D(2)));
                res=max(res,abs(dC2/Ctau(2)));
                res=max(res,abs(dU2/U(2)));
                if res>0.3; Rel=0.3/res; else Rel=1; end
                T2=T(2) + Rel*dT2;
                D2=D(2) + Rel*dD2;
                C2=Ctau(2) + Rel*dC2;   
                U2=U(2) + Rel*dU2;
                U(2)=U2;
            end
            %----------------------------------------------------------------------------------
            
            % filter extrem values
            C2=min(C2,0.3);
            C2=max(C2,0.0000001);
            
            if IsWake; D2=D2 + LTE(i+step); Hmin=1.00005;else Hmin=1.02; end
            % filter to small H
            dh= max(0,Hmin-D2/T2);
            D2=D2 +dh*T2;
            
            k=k+1; % iteration counter
        end % <- End of iteration loop for current intervall
        % update Ctau solution
        sol.c(i+step)=C2;C1=C2;
        % update parameter values
        sol.UQ(i)=params(1,7);
        sol.Us(i)=params(1,8);
        sol.CtEQ(i)=params(1,9);
        sol.Del(i)=params(1,10);
        
    end 
    
%------------------------------------------------- Update panel point 2 solution -------------------------------------       
    % Filter if solution does not make sense -> extrapolate solution
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
            T2=T1*sqrt(s(i+step)/s(i));
            D2=D1*sqrt(s(i+step)/s(i)); 
        end
        if i~=Ende; U2= (sol.U(i+step)+sol.U(i-step))/2;end
        sol.HK(i+step)=D2/T2;
    end
    %---------------------------------------------
    
    % store calculated 2 value
    sol.T(i+step)=T2;
    sol.D(i+step)=D2;
    sol.U(i+step)=U2;

    % store parameters
    sol.CfM(i)=params(1,1);
    sol.Cf(i)=params(1,2);
    sol.CDM(i)=params(1,3);
    sol.CD(i)=params(1,4);
    sol.HS(i)=params(1,5);
    sol.Ret(i)=params(1,6);
    
end  % <- End of the loop over all boundary subintervalls

% last values
   sol.Cf(i+step)=params(2,2);
   sol.CD(i+step)=params(2,4);
   sol.HS(i+step)=params(2,5);
   sol.Ret(i+step)=params(2,6);
   if length(params(1,:))>7 % additional turbulent parameters
       sol.UQ(i+step)=params(2,7);
       sol.Us(i+step)=params(2,8);
       sol.CtEQ(i+step)=params(2,9);
       sol.Del(i+step)=params(2,10);
   end
   
end






