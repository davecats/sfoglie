function [ TM,T2,DM,D2,U2,UM,C2,Llam] = Transition(sol,i,step,Cst,U,Vb,h)
%TRANSITION finds transition point on Intervall of transition. Splits Intervall in laminar and turbulent part
%           and calculates the values of delta1 and delta2 for the transition point and the panels 2 point

nu=evalin('base','nu');
% load quantities
n1=sol.c(i);
n2=sol.c(i+step);
T1=sol.T(i);
D1=sol.D(i);
HK1=sol.HK(i);
HKset=false;

% approximates  transition point with secant method
Llam= (9-n1)/(n2-n1)*h;
Lturb= h-Llam;

% velocities at transitionpoint
UM=(U(2)*Llam + U(1)*Lturb)/h;
VM=(Vb(2)*Llam + Vb(1)*Lturb)/h;

%---------------------------------------- laminar part----------------------------------------------------
% calculate values at transition point with laminar treatment and gets exact transition length

k=0; res=1;
D2=D1;T2=T1; 
Ul=[U(1);UM];
Vl=[Vb(1);VM];
while res>sol.resmax && k<sol.itmax
    Dl=[D1;D2];Tl=[T1;T2];

    [ f1,f2,df_dT,df_dD, ~,~,df_dh ] = SingleJacobiLam(Dl,Tl,Ul, Vl,Llam,[0;0],true,false); 
    % add line for Amplification factor that forces n=ncrit and adds laminar length as variable

    [dn ,dn_dT,dn_dH,dn_dRet] = AmplificationDerivate(Dl./Tl,Tl.*Ul/nu,Tl );
    dn_dT2=dn_dT(:,2) - dn_dH(:,2)*Dl(2)/Tl(2)^2 + dn_dRet(:,2)*Ul(2)/nu;
    dn_dD2=           + dn_dH(:,2)*1/Tl(2);
    J=[df_dT, df_dD, df_dh; -Llam*dn_dT2, -Llam*dn_dD2, -dn];
    tmp=n1+Llam*dn-9; % forces nT=9
    rhs=[-f1;-f2; tmp];
    dz=J\rhs;
    dT2=dz(1);
    dD2=dz(2);
    dL =dz(2);

    res=max(abs(dT2/T2),abs(dD2/D2));
    res=max(res,abs(dL/Llam));   
    % under relaxation for big changes
    if res>0.3; Rel=0.3/res; else Rel=1; end
    % update values
    T2=T2+Rel*dT2;
    D2=D2+Rel*dD2;
    Ltmp=Llam+Rel*dL;
    
    HK2=D2/T2;
    % seperation when shape parameter reached a certain value
    % -> prescribe th growth of H and adapts the tangential boundary edge velocity
    %----------------------------------------------------------------------------------

    if HK2 > sol.HmaxLam 
        HKset=true;
        Htmp=HK1 + 0.03*Llam/Tl(1) ; % limit increase of H
        Htmp=max(Htmp,sol.HmaxLam);

        %if k==0;disp(['HK set at transition point HK=' num2str(Htmp)  ] );end

        Dl(2)=Htmp*Tl(2); % correct displacement thickness

        % add condition H2 = Htmp to Newton System and solves for new velocity U2
        [ f1,f2,df_dT,df_dD,~,df_dU,df_dh  ] = SingleJacobiLam(Dl,Tl,Ul, Vl,Llam,[0; 0],true,false); 
        
        [dn ,dn_dT,dn_dH,dn_dRet] = AmplificationDerivate(Dl./Tl,Tl.*Ul/nu,Tl);
        dn_dT2=dn_dT(:,2) - dn_dH(:,2)*Dl(2)/Tl(2)^2 + dn_dRet(:,2)*Ul(2)/nu;
        dn_dD2=           + dn_dH(:,2)*1/Tl(2);
        dn_dU2=                                      + dn_dRet(:,2)*Tl(2)/nu;
        
        J=[df_dT, df_dD, df_dU,df_dh; -Dl(2)/Tl(2)^2, 1/Tl(2),0,0 ; -Llam*dn_dT2, -Llam*dn_dD2, -Llam*dn_dU2, -dn];
        tmp=n1+Llam*dn-9;
        rhs=[-f1;-f2; Htmp - Dl(2)/Tl(2) ;tmp];
        dz=J\rhs;
        dT2=dz(1);dD2=dz(2);dU2=dz(3);dL=dz(4);
        res=max(abs(dT2/Tl(2)),abs(dD2/Dl(2)));
        res=max(res,abs(dU2/Ul(2)));
        res=max(res,abs(dL/Llam));
        if res>0.3; Rel=0.3/res; else Rel=1; end
        T2=Tl(2)+Rel*dT2;
        D2=Dl(2)+Rel*dD2;
        U2=Ul(2)+Rel*dU2;
        Ltmp=Llam+Rel*dL;
        Ul(2)=U2;               
    end
    %----------------------------------------------------------------------------------            
    Llam=Ltmp;
    dh= max(0,1.02-D2/T2);
    D2=D2 +dh*T2;

    % update transition point value for V
    VM=(Vb(2)*Llam + Vb(1)*Lturb)/h;
    Vl(2)=VM;
    
    k=k+1;
%     if k==sol.itmax;
%         disp(['not konverged on transition panel, residuum: ' num2str(res)] );
%     end
end

% test
% a=IntegrateAmplification(n1,Llam,T1,T2,D1/T1,D2/T2,Ul(1)*T1/nu,U(2)*T2/nu); 

% save new laminar length
Lturb=h-Llam;


C2=Cst;
% calculate start value for Ctau from equilibrium shear stress
fac=1.8*exp(-3.3./(HK1-1));
%fac=1.1*exp(-10./HK2.^2);
[ HS, ~, ~  ]=H32turb( Dl./Tl,Ul.*Tl/nu);
Us=0.5*HS.*( -1/3 + Tl./(0.75*Dl) );
[CEQ,~,~,~,~]=CtEQ( Dl./Tl,Ul.*Tl/nu,HS,Us,false);
C1=fac* CEQ(1);


if HKset
    UM=Ul(2); 
    HK1= HK2; 
end


TM=T2; DM=D2;
D1=DM; T1=TM;

%---------------------------------- turbulent part----------------------------------------------------


% calculate values at 2 point with turbulent treatment
k=0; res=1;
D2=D1;T2=T1; 
Ut=[UM;U(2)];
Vt=[VM,Vb(2)];


Dt=[D1;D2];Tt=[T1;T2];

 while res>sol.resmax && k< sol.itmax
        Ctau=[C1;C2];

        [ f1,f2,f3,df_dT,df_dD,df_Ct,~ ] = SingleJacobiTurb(Dt,Tt,Ctau,Ut,Vt,Lturb,[0; 0],false,false); 

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
        % seperation when shape parameter reached a certain value
        % -> prescribe th growth of H and adapts the tangential boundary edge velocity
        %----------------------------------------------------------------------------------
        if HK2>sol.HmaxTurb   
            Htmp= HK1-0.15*Lturb/T1;
            Htmp=max(Htmp,sol.HmaxTurb);
            %if k==0;disp(['HK set at node ' num2str(i+step) ' HK=' num2str(Htmp)  ] );end

            Dt(2)=Htmp*Tt(2); % correct displacement thickness

            % add condition H2 = Htmp to Newton System and solves for new velocity U2
            [ f1,f2,f3,df_dT,df_dD,df_Ct,~,df_dU] = SingleJacobiTurb(Dt,Tt,Ctau,Ut,Vb,Lturb,[0; 0],true,false);
            J=[df_dT, df_dD, df_Ct,df_dU; 0, -Dt(2)/Tt(2)^2, 1/Tt(2),0];
            rhs=[-f1;-f2;-f3; Htmp- Dt(2)/Tt(2)];
            dz=J\rhs;
            dT2=dz(1);dD2=dz(2);dC2=dz(3);dU2=dz(4);
            res=max(abs(dT2/Tt(2)),abs(dD2/Dt(2)));
            res=max(res,abs(dC2/Ctau(2)));
            res=max(res,abs(dU2/Ut(2)));
            if res>0.3; Rel=0.3/res; else Rel=1; end
            T2=Tt(2)+Rel*dT2;
            D2=Dt(2)+Rel*dD2;             
            U2=Ut(2)+Rel*dU2;
            C2=Ctau(2)+Rel*dC2;
            Ut(2)=U2;
        end
        Tt(2)=T2; Dt(2)=D2;
        
        k=k+1;
        %if k==sol.itmax;disp(['not konverged on transition panel, residuum: ' num2str(res)] );end
 end
        

end

