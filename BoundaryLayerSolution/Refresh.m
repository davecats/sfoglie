function sol= Refresh(prf,flo,sol,eng)
%Refresh goes through the Boundary layer after each Newton iteration to
%        find the new transition point and force a certain H-U-correlation
%        for the seperation region
     

nkrit=flo.nkrit;
nu=flo.nu;

Lges=[prf.panels.L, flo.wake.L]; 

% initial node values
[ sol.T(prf.Nle-1), sol.D(prf.Nle-1) ] = ImproveStartNodeSol( nu,sol.U(prf.Nle-1),prf.LE1,sol.Vb(prf.Nle-1));
[ sol.T(prf.Nle)  , sol.D(prf.Nle)   ] = ImproveStartNodeSol( nu,sol.U(prf.Nle)  ,prf.LE2,sol.Vb(prf.Nle)  );

% initial node values
sol.c(prf.Nle-1)=0;
sol.c(prf.Nle)=0;

sol.HK(prf.Nle-1) =sol.D(prf.Nle-1)/sol.T(prf.Nle-1);
sol.Ret(prf.Nle-1)=sol.T(prf.Nle-1)*sol.U(prf.Nle-1)/nu;

sol.HK(prf.Nle) =sol.D(prf.Nle)/sol.T(prf.Nle);
sol.Ret(prf.Nle)=sol.T(prf.Nle)*sol.U(prf.Nle)/nu;


for section=1:3
    if section== 1 % suction side
        lam=true; Isflo.wake=false;
        Start=prf.Nle-1; step=-1; Ende=2; shift=1;
        s=prf.sU; 
        sBL=s;% arc length in boundary layer direction
        gap=zeros(prf.Nle-1,1);
        TranOld=Start - sol.iTran(1);
        j=1; % BL node index
        tr=false;
    elseif section == 2  % pressure side
        lam=true; Isflo.wake=false;
        Start=prf.Nle; step=1; Ende=prf.N-1; shift=0;
        s=prf.s; 
        sBL=[zeros(1,prf.Nle-1),prf.sL];
        gap=zeros(prf.N,1);
        TranOld=sol.iTran(2) - Start;
        j=1;
        tr=false;
    else  % wake
       lam=false;  Isflo.wake=true;
       Start=prf.N +1; step=1; Ende=prf.N+flo.wake.N-1; shift=0;
       s=flo.wake.s; 
       sBL= [ zeros(1,prf.N),flo.wake.s(1:end-1)+prf.sL(end)]; 
       gap=[zeros(prf.N,1);flo.wake.gap];
       j=prf.N+1;
       sol.T(j)=sol.T(1)+sol.T(prf.N);
       sol.D(j)=sol.D(1) + sol.D(prf.N) + prf.gap; 
       sol.HK(j)=(sol.D(j)-prf.gap)/sol.T(j); % TE gap does not go into shape parameter
    end
    
    
    %--------------------------------------------------------------------
    
    for i=Start:step:Ende
        ind=i:step:i+step;
        Tran=max(Start,Ende);
        
        
        H2=sol.HK(ind(2));
        U2=sol.U(ind(2));
        T2=sol.T(ind(2));
        D2=sol.D(ind(2));
        C2=0.03;
        
        % values of direct mode as reference values
        Uref=U2;
        Href=H2;
        
        if j>= TranOld % former turbulent node is now laminar
            C2=sol.c(i+step); 
            if C2<0; C2=0.03; end  
            if lam
                Href=( sol.D(i)-gap(i) )/sol.T(i); % Take value of previous node as reference
            end  
        elseif j>=Tran && j<TranOld % former laminar node is now turbulent
            if j==Tran
                C2=0.03;%InitialCtau(D2,T2,U2,nu);
            else
                C2=sol.c(i); % Take value of previous node as reference -> already known
            end
            if C2<0; C2=0.03; end

        end

        
        res=1; k=0;
        % ------ iteration loop -------
        while res>5e-6 && k<sol.itmax
            D=[sol.D(i);D2];
            T=[sol.T(i);T2];
            U=[sol.U(i);U2];
            H= (D-gap(ind))./T;
            Ret=U.*T/nu;
            
            if lam     
                [ f1,f2,df_dT,df_dD,df_dU ] = SingleJacobiLam(T,U, sol.Vb(ind),H,Ret,Lges(i-shift),sBL(i),true,nu,sol.pressureTerm(ind)); 
                
                % calculate changes in U for the bigger shape parameter H + 1 where the equation is still fullfilled
                J=[df_dT, df_dD, df_dU; -H2./T2, 1./T2 ,0];
                rhs=[-f1;-f2;1];
                
                d=J\rhs;
                if k<16 % sensitivity of U to the change of H
                    % constant factor determines the allowed changes in H
                    %   -> bigger values lead to bigger allowed deviation of H
                    senTMP= 1000 * d(3)* Href/Uref;
                    if k<6; sen=senTMP; else sen= 0.5*(sen+senTMP); end
                end
                
                % calculate variables for the forced change in H
                J(3,1)=J(3,1)*Href;
                J(3,2)=J(3,2)*Href; 
                J(3,3)=sen/Uref;
                rhs(3)=-Href^2*(H2/Href -1) - sen*(U2/Uref -1);
                
                dz=J\rhs;
                
                if ~isempty(find(~isreal(dz),1)); res=1; break; end % prevent complex values
                
                res=max(abs( [dz(1)/T2,dz(2)/D2,dz(3)/U2] ));
                % under relaxation for big changes
                if res>0.3; Rel=0.3/res; else Rel=1; end            
                % update values
                T2=T2 + Rel*dz(1);
                D2=D2 + Rel*dz(2);
                U2=U2 + Rel*dz(3);

            else % ---------------- turbulent -----------------
                
                [ f1,f2,f3,df_dT,df_dD,df_Ct,df_dU] = ...
                    SingleJacobiTurb(D,T,[sol.c(i);C2],U, sol.Vb(ind),H,Ret,Lges(i-shift),sBL(i),true,Isflo.wake,gap(ind),nu,sol.pressureTerm(ind));

                J=[df_dT, df_dD, df_Ct,df_dU; -H2./T2, 1./T2,0,0];
                rhs=[-f1;-f2;-f3; 1];
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
                
                if ~isempty(find(~isreal(dz),1)); res=1; break; end % prevent complex values
                
                res=max(abs( [dz(1)/T2,dz(2)/D2,dz(3)/C2,dz(4)/U2] ));

                if res>0.3; Rel=0.3/res; else Rel=1; end
                T2=T2 + Rel*dz(1);
                D2=D2 + Rel*dz(2);
                C2=C2 + Rel*dz(3);   
                U2=U2 + Rel*dz(4);  
                
                % filter extrem values
                C2=min(C2,0.3);
                C2=max(C2,0.0000001);
            end
             
            % correction for to small H values
            if Isflo.wake; Hlim=1.00005; else Hlim=1.02; end
            dh= max(0,Hlim - ( D2-gap(i+step) )/T2);
            D2= D2 +dh*T2;
            H2= ( D2-gap(i+step) )/T2;
           
           k=k+1; 
        end
        sen=senTMP;
        
        if res>0.1 
            %disp(['not konverged->solution extrapolated, node: ' num2str(i+step) ', residuum: ' num2str(res)] );
            if Isflo.wake
                T1=sol.T(i); D1=sol.D(i);
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
                if i>Tran; C2= sol.c(i); end
            end
            if  i~=Ende; % central approximation
                U2= (sol.U(i+2*step)*Lges(i-shift)+sol.U(i)*Lges(i-shift+step))/(Lges(i-shift+step)+Lges(i-shift));
            else % backward approximation
                U2= sol.U(i);
            end
            sol.HK(i+step)=( D2-gap(i+step) )/T2;
        end
        
        if T2<0 || D2<0 % in case of negativ values -> do not refresh
             T2=sol.T(i+step);
             D2=sol.D(i+step);
             C2=sol.c(i+step);
             U2=sol.U(i+step);
        end

        if lam % ----------- Transition
            
            n2 = AmplSol(flo,sol.c(i),[sol.T(i);T2],[sol.U(i);U2], [sol.HK(i); D2/T2],[sol.Ret(i); T2*U2/nu], Lges(i-shift) );
            sol.c(i+step)=n2;
            
            
             % tripping arc length in current intervall
            if sol.Tripping(section) && prf.nodes.X(i+step)>sol.xT(section) && prf.nodes.X(i)<sol.xT(section);
                tr=true;
                w1= (prf.nodes.X(i+step)-sol.xT(section))/ (prf.nodes.X(i+step)-prf.nodes.X(i));
                w2=1-w1;
                sol.sT(section)= w1*sBL(i) + w2*sBL(i+step);
            end


            % transition occures in current intervall
            if sol.c(i+step)>nkrit || tr; 
                
                % free transition before tripping
                if ~tr && sol.Tripping(section)
                    sol.Tripping(section)=false; % set tripping false for transition EQ
                    before=true;
                else
                    before=false;
                end 
                
                lam=false;  
                % solve the Equations for Transition panel 
                [ tran,T2,D2,U2,C1,~ ] = RefreshTransition(section, flo,eng,sol,sol.c(ind), sol.T(ind), sol.D(ind),sol.U(ind),sol.Vb(ind),sBL(i),Lges(i-shift),C2 );
               
                
                if before; sol.Tripping(section)=true; end % make sure tripping location is still checked
                
                %----------
                % write transition node values
                sol.iTran(section)=i+step;
                sol.tran.s(section)=prf.s(i)-tran.Llam;
                sol.tran.n2(section)=sol.c(i+step);
                sol.tran.Llam(section)=tran.Llam;
                sol.tran.Lturb(section)=tran.Lturb;
                sol.tran.wl(section)=tran.wl;
                sol.tran.wt(section)=tran.wt;
                sol.tran.x(section)=tran.wl*prf.nodes.X(i)+tran.wt*prf.nodes.X(i+step);
                
                % store Ctau instead of n 
                C2=C1;
                sol.c(i+step)=C1; 
            end

            %if no transition occured so far -> set transition to last airfoil node
            if i==Ende
                lam=false;
                sol.Tripping(section)=true;
                sol.sT(section)=sBL(i+step);
                sol.xT(section)=prf.nodes.X(i+step);
                sol.tran.x(section)=prf.nodes.X(i+step);
                
                % set initial value for sqrt(Ctau) at TE
                C= InitialCtau( sol.D(ind),sol.T(ind),sol.U(ind),nu );
                C1=C(2);
                %----------            
                % write transition node values
                sol.iTran(section)=i+step;
                sol.tran.s(section)=prf.s(i+step);
                sol.tran.x(section)=prf.nodes.X(i+step);
                sol.tran.n2(section)=sol.c(i+step);
                sol.tran.Llam(section)=Lges(i-shift);
                sol.tran.Lturb(section)=0;
                
                C2=C1;
                sol.c(i+step)=C1; 
                T2=sol.T(i+step);D2=sol.D(i+step);U2=sol.U(i+step);
                
            end 
            %------------------------------------------------------------------------------------------
        end
        
         % Write new values if they are plausible  
         if T2>0 && D2>0 % in case of negativ values -> do not refresh
            if ~lam; sol.c(i+step)=C2; end 
            sol.T(i+step)=T2;
            sol.D(i+step)=D2;
            sol.U(i+step)=U2;
            sol.HK(i+step)= (D2-gap(ind(2))) /T2;
            if sol.HK(i+step)<0
               disp(['neg. D', num2str(D2), num2str(T2)]) 
            end
            
            
            sol.Ret(i+step)=U2*T2/nu;
         end
        
        j=j+1;
    end 
    
    
    
end

if prf.sharpTE
   % make sure the velocities for both TE points are the same 
   sol.U(1)=sol.U(prf.N+1);
   sol.U(prf.N)=sol.U(prf.N+1);
end


end






