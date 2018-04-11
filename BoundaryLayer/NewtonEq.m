function [solnew, prf] = NewtonEq( prf,wake,sol,D0,Uinv,it)
%NEWTONEQ       sets up the global Newton equationsystem J dz=- f(z) and
%               solves it iterativly
%               sol:    initial solution struct
%               D0:      Coefficients of mass Defekt -> U = Uinv + Dm
%               it:     max number of iteration steps

%nu=evalin('base','nu');
%Re=evalin('base','Re');

Nges=prf.N + wake.N;
k=0;
res=1;


% derivates of stagnation point position in respect to Gamma
prf.dsLE_dG1= prf.LE2/(Uinv(prf.Nle) + Uinv(prf.Nle-1));
prf.dsLE_dG2= prf.LE1/(Uinv(prf.Nle) + Uinv(prf.Nle-1));

% save start leading edge position
%NleStart=prf.Nle;
GamStart= [Uinv(1:prf.Nle-1); -Uinv(prf.Nle:end)];

while res>3e-4 && k<it 
     % correct adjust sign for pressure/suction side
     sgn=ones(size(D0));
     sgn(1:prf.Nle-1,1:prf.Nle-1)=-sgn(1:prf.Nle-1,1:prf.Nle-1);
     sgn(prf.Nle:end,prf.Nle:end)=-sgn(prf.Nle:end,prf.Nle:end);
     
     D=D0.*sgn;   
    
     if k>0;
        sol=solnew;
        % refresh -> search new transition point etc
        sol = Refresh( prf,wake, sol );
     end
     Un=Uinv + D*sol.m; % velocity solution with new source Terms

    [J,rhs]=JacobiM(prf,wake,sol, Un,D);
    dz=J\rhs; % Solution of equation -> gives correction for current Newton-step
    
    dT=dz(1:Nges);
    dc=dz(Nges+1:2*Nges);
    dm=dz(2*Nges+1:end);

    % updates values with underrelaxation for big changes
    [solnew,Rel,res]=Update(prf,wake,sol,dT,dc,dm,Uinv,D,k);

    % max residuum
    k=k+1;
    
    disp(['Iteration step ' num2str(k) ' Relaxation factor ' num2str(Rel) ' max residual ' num2str(res) ] );
    
    % Get new stagnation point and adjust values
    gam=[solnew.U(1:prf.Nle-1);  -solnew.U(prf.Nle:end)];
    NleOld=prf.Nle;
    % find new stagnation point
    [ prf.Nle,prf.sLE,prf.LE1,prf.LE2 ] = getStagnationPoint( gam, prf.s );
    prf.sL=prf.s(prf.Nle:end)-(prf.s(prf.Nle)-prf.LE2)*ones(1,prf.N-prf.Nle+1); 
    prf.sU=(prf.s(prf.Nle-1)+prf.LE1)*ones(1,prf.Nle-1)-prf.s(1:prf.Nle-1); 

    % stagnation point position derivates in respect to gamma
    prf.dsLE_dG1= prf.LE2/(solnew.U(prf.Nle) + solnew.U(prf.Nle-1) );
    prf.dsLE_dG2= prf.LE1/(solnew.U(prf.Nle) + solnew.U(prf.Nle-1) );
    
    NleDif=prf.Nle-NleOld; 
    % set new startinpoint variables
    if NleDif>0 %-> moves more on pressure side
        solnew.T(NleOld:prf.Nle-1)=solnew.T(prf.Nle);
        solnew.D(NleOld:prf.Nle-1)=solnew.D(prf.Nle);
        solnew.c(NleOld:prf.Nle-1)=solnew.c(prf.Nle);
        UL=solnew.U(NleOld-1)/prf.sU(end-NleDif);
        st=prf.sU(end:-1:1);
        solnew.U(NleOld:prf.Nle-1)= UL*st(1:abs(NleDif));        
        disp(['Stagnationpoint moved: Old ',num2str(NleOld),' new ',num2str(prf.Nle)])
        % adjust inviscid velocity
        Uinv=sign(gam).*GamStart;
    elseif NleDif<0 %-> moves more on suction side
        solnew.T(prf.Nle:NleOld-1)=solnew.T(NleOld);
        solnew.D(prf.Nle:NleOld-1)=solnew.D(NleOld);
        solnew.c(prf.Nle:NleOld-1)=solnew.c(NleOld);
        UL=solnew.U(NleOld)/prf.sL(1-NleDif);
        solnew.U(prf.Nle:NleOld-1)= UL*prf.sL(1:abs(NleDif));
        disp(['Stagnationpoint moved: Old ',num2str(NleOld),' new ',num2str(prf.Nle)])
        % adjust inviscid velocity
        Uinv=sign(gam).*GamStart;
    end
    % make shure the velocitie does not get zero
    solnew.U(solnew.U<1e-7)=1e-7;
    solnew.m=solnew.D.*solnew.U;    
    
    % 
    
end

%solnew = Refresh( prf,wake, solnew );

% final values
indL=(solnew.iTran(1)+1:solnew.iTran(2)-1); % laminar node indizes
indT=[1:solnew.iTran(1) , solnew.iTran(2):prf.N]; % wake node indizes
indW=prf.N+1:prf.N+wake.N; % turbulent node indizes

I=zeros(size(solnew.T));
solnew.Cf=I;
solnew.Cf(indL)= CF2lam(solnew.HK(indL), solnew.Ret(indL));
solnew.Cf(indT)= CF2turb(solnew.HK(indT), solnew.Ret(indT));

solnew.tau=solnew.Cf .*(solnew.U/prf.Uinfty).^2;
solnew.Cf=2*solnew.Cf;

solnew.HS=I;
solnew.HS(indL)=H32lam(solnew.HK(indL));
solnew.HS(indT)=H32turb(solnew.HK(indT),solnew.Ret(indT));
solnew.HS(indW)=H32turb(solnew.HK(indW),solnew.Ret(indW));

solnew.Us=0.5*solnew.HS.*( -1/3 + 1./(0.75*solnew.HK) );


solnew.CD=I;
solnew.CD(indL)=CD2lam(solnew.HK(indL), solnew.Ret(indL));
solnew.CD(indT)=CD2turb(solnew.HK(indT),solnew.Ret(indT),solnew.HS(indT), solnew.Us(indT),solnew.c(indT),false);
solnew.CD(indW)=CD2turb(solnew.HK(indW),solnew.Ret(indW),solnew.HS(indW), solnew.Us(indW),solnew.c(indW),true);

solnew.CD=0.5*solnew.HS.*solnew.CD;


solnew.Cp=1-(solnew.U/prf.Uinfty).^2;
if prf.IsSharp
   solnew.Cp(prf.N)=solnew.Cp(1); 
end



Xt=prf.nodes.X(solnew.iTran);

tmp= solnew.iTran +[1,-1];
Xl=prf.nodes.X(tmp);

wl= solnew.tran.Lturb./(solnew.tran.Llam + sol.tran.Lturb)  ;
wt= solnew.tran.Llam./(solnew.tran.Llam + sol.tran.Lturb)  ;

xTran= wl.*Xl + wt.*Xt;

% Write out
if ~solnew.Tripping(1); 
    solnew.xT(1)=xTran(1); 
    disp(['suction side: transition at x/c=',num2str(xTran(1))]);
else
    disp(['suction side: forced transition at x/c=',num2str(solnew.xT(1))]);
end
if ~solnew.Tripping(2); 
    solnew.xT(2)=xTran(2); 
    disp(['pressure side: transition at x/c=',num2str(xTran(2))]);
else
    disp(['pressure side: forced transition at x/c=',num2str(solnew.xT(2))]);
end  

% Integral values

[solnew.CL,solnew.Cdrag,solnew.Cnu]=IntValues(solnew.Cp(1:prf.N),solnew.tau(1:prf.N),prf.s,prf.panels.e',prf.panels.n',prf.alfa,...
    true,solnew.Vb/prf.Uinfty);


%old
% % lift coefficient
% h=( prf.panels.X(2,:)-prf.panels.X(1,:) )*cos(prf.alfa) + ( prf.panels.Y(2,:)-prf.panels.Y(1,:) )*sin(prf.alfa);
% h=h';
% 
% Cpp=[solnew.Cp(1:prf.N);solnew.Cp(1)];
% solnew.CL   = 0.5*sum( h.*( Cpp(2:end) + Cpp(1:end-1) ) );% midpoint rule
% 
%     % shear forces of each panel in x and y direction 
% tauX= abs(solnew.tau(1:prf.N).*prf.panels.e(1,:)');
% tauY= abs(solnew.tau(1:prf.N).*prf.panels.e(2,:)');
% 
% tauXU= [0;tauX(prf.Nle-1:-1:1)];
% tauYU= [0;tauY(prf.Nle-1:-1:1)];
% tauXL= [0;tauX(prf.Nle:prf.N)];
% tauYL= [0;tauY(prf.Nle:prf.N)];
% sU=[0,prf.sU(end:-1:1)];
% sL=[0,prf.sL];
% 
% % integrate x-component -> simpson law
% FX=NumInt(tauXU,sU) + NumInt(tauXL,sL);
% FY=NumInt(tauYU,sU) + NumInt(tauYL,sL);
% 
% % viscous Drag Coefficient
% solnew.Cnu=2*( FX*cos(prf.alfa) + FY*sin(prf.alfa) );
% 
% % Drag coefficient 
% solnew.Cdrag=DragCoeff(solnew.T(end),solnew.HK(end),solnew.U(end)/prf.Uinfty, 1 );


% integral Cf
indU=prf.Nle-1:-1:1;
indL=prf.Nle:prf.N;

% midpoint rule
[~,tmp]=NumInt([0;solnew.tau(indU)],[0,prf.sU(end:-1:1)] ,2,true);
solnew.CI_U=tmp;
[~,tmp]=NumInt([0;solnew.tau(indL)],[0,prf.sL]           ,2,true);
solnew.CI_L=tmp;


end


