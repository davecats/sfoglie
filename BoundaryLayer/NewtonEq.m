function [solnew, prf] = NewtonEq( prf,flo,eng,sol,D0,Uinv,it,pressureCor)
%NEWTONEQ       sets up the global Newton equationsystem J dz=- f(z) and
%               solves it iterativly
%               sol:    initial solution struct
%               D0:     Coefficients of mass Defekt -> U = Uinv + Dm
%               it:     max number of iteration steps



Nges=prf.N + flo.wake.N;
k=0;
res=1;


% derivates of stagnation point position in respect to Gamma
prf.dsLE_dG1= prf.LE2/(Uinv(prf.Nle) + Uinv(prf.Nle-1));
prf.dsLE_dG2= prf.LE1/(Uinv(prf.Nle) + Uinv(prf.Nle-1));

% save start leading edge position
%NleStart=prf.Nle;
GamStart= [Uinv(1:prf.Nle-1); -Uinv(prf.Nle:end)];



% Model for correction of the pressure Term   
%------------------------------------------------------------------------------------------
if nargin==7
    pressureCor=false;
end
Blow=find(sol.Vb~=0);

if ~isempty(Blow) && pressureCor
     BlowU= Blow(Blow<prf.Nle-1);
     BlowL= Blow(Blow>prf.Nle);
     
     pressureTerm=zeros(size(sol.Vb));
     if ~isempty(BlowU)
        sU_b=prf.sU(BlowU);
        PrU = PressureCorrect(sU_b(end),sU_b(1),prf.sU(end:-1:1)',sol.Vb(prf.Nle-1:-1:1),sol.U(prf.Nle-1:-1:1));
        PrU=PrU(end:-1:1);
        pressureTerm(1:prf.Nle-1)=PrU;
     end
     if ~isempty(BlowL)
        sL_b=prf.sL(BlowL-prf.Nle+1);
        PrL = PressureCorrect(sL_b(1),sL_b(end),prf.sL',sol.Vb(prf.Nle:prf.N),sol.U(prf.Nle:prf.N));
        pressureTerm(prf.Nle:prf.N)=PrL;
     end
     % add effect of pressure Term to Vb     
     Vb_base=sol.Vb;
     sol.Vb=sol.Vb - pressureTerm;
end
%------------------------------------------------------------------------------------------



%xtrr=zeros(it,2);

while res>eng.tol && k<it 
    
     % adjust sign for pressure/suction side
     sgn=ones(size(D0));
     sgn(1:prf.Nle-1,1:prf.Nle-1)=-sgn(1:prf.Nle-1,1:prf.Nle-1);
     sgn(prf.Nle:end,prf.Nle:end)=-sgn(prf.Nle:end,prf.Nle:end);
     
     D=D0.*sgn;   
    
     if k>0;
        sol=solnew;
        % refresh -> search new transition point and runs inverse mode for seperation
        sol = Refresh( prf,flo, sol, eng );
     end
     %xtrr(k+1,:)=sol.tran.x';
     
     Un=Uinv + D*sol.m; % velocity solution with new source Terms

    [J,rhs]=JacobiM(prf,flo,eng,sol, Un,D);
    dz=J\rhs; % Solution of equation -> gives correction for current Newton-step
    
    dT=dz(1:Nges);
    dc=dz(Nges+1:2*Nges);
    dm=dz(2*Nges+1:end);

    % updates values with underrelaxation for big changes
    [solnew,Rel,res]=Update(prf,flo,sol,dT,dc,dm,Uinv,D,k);

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
    % set new starting point variables
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
    % make shure the velocity does not get zero
    solnew.U(solnew.U<1e-7)=1e-7;
    solnew.m=solnew.D.*solnew.U;    
    
    % 
    
end

%solnew = Refresh( prf,flo.wake, solnew );

solnew.residual=res;


% figure
% hold on
% plot(1:it,xtrr(:,1),'b .')
% plot(1:it,xtrr(:,2),'r x')
% xlabel('iter')
% ylabel('x_tran')
% legend('suction','pressure')

if ~isempty(Blow) && pressureCor
    % substract effect of pressure to get Normal vb Term   
     solnew.Vb=Vb_base;
     solnew.pressureTerm=pressureTerm;
end


% final values
% calculate values of end solution (Cf, CD, CL ...)

indL=(solnew.iTran(1)+1:solnew.iTran(2)-1); % laminar node indizes
indT=[1:solnew.iTran(1) , solnew.iTran(2):prf.N]; % flo.wake node indizes
indW=prf.N+1:prf.N+flo.wake.N; % turbulent node indizes

I=zeros(size(solnew.T));
solnew.Cf=I;
solnew.Cf(indL)= CF2lam(solnew.HK(indL), solnew.Ret(indL));
solnew.Cf(indT)= CF2turb(solnew.HK(indT), solnew.Ret(indT));

solnew.tau=solnew.Cf .*(solnew.U/flo.Uinfty).^2;
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


solnew.Cp=1-(solnew.U/flo.Uinfty).^2;
if prf.sharpTE
   solnew.Cp(prf.N)=solnew.Cp(1); 
end

% find seperation points
% suctionside
[solnew.xseparation(1), solnew.xreattach(1)]=FindSeparationLoc( prf.nodes.X(prf.Nle-1:-1:1),solnew.HK(prf.Nle-1:-1:1),...
                                                                                            prf.Nle-solnew.iTran(1), 3.8,2.5 );
%pressure side
[solnew.xseparation(2), solnew.xreattach(2)]=FindSeparationLoc( prf.nodes.X(prf.Nle:prf.N),solnew.HK(prf.Nle:prf.N),...
                                                                                           solnew.iTran(2)-prf.Nle+1, 3.8,2.5 );


% Write out
disp(' ')
if solnew.Tripping(1); 
    disp(['suction side: forced transition at x/c=',num2str(solnew.xT(1))]);
else
    disp(['suction side: free transition at x/c=',num2str(solnew.tran.x(1))]);
end
if solnew.Tripping(2); 
    disp(['suction side: forced transition at x/c=',num2str(solnew.xT(2))]);
else
    disp(['suction side: free transition at x/c=',num2str(solnew.tran.x(2))]);
end
disp(' ')

% Integral values

 if prf.sharpTE
    [solnew.CL,solnew.Cdrag,solnew.Cnu]=IntValues(solnew.Cp(1:prf.N-1),solnew.tau(1:prf.N-1),prf.s(1:prf.N-1),prf.nodes.e(:,1:prf.N-1)',...
                                         prf.nodes.n(:,1:prf.N-1)', flo.alfa, true,solnew.Vb/flo.Uinfty);
 else
%     [solnew.CL,solnew.Cdrag,solnew.Cnu]=IntValues(solnew.Cp(1:prf.N),solnew.tau(1:prf.N),prf.s,...
%                                                 prf.nodes.e',prf.nodes.n',flo.alfa, true,solnew.Vb/flo.Uinfty);
    %include TE Panel contribution
    [solnew.CL,solnew.Cdrag,solnew.Cnu]=IntValues(solnew.Cp(1:prf.N+1),solnew.tau(1:prf.N+1),[prf.s,prf.s(end)+prf.panels.L(end)],...
                                                [prf.nodes.e,prf.panels.e(:,end)]',[prf.nodes.n,prf.panels.n(:,end)]',...
                                                flo.alfa, true,solnew.Vb/flo.Uinfty);
 end

if prf.N<80 % numeric integration for drag gets to inaccurate -> use Squire-Young formula
    solnew.Cdrag=DragCoeff(solnew.T(end),solnew.HK(end),solnew.U(end)/flo.Uinfty, 1 );
end

% integral Cf
indU=prf.Nle-1:-1:1;
indL=prf.Nle:prf.N;

% midpoint rule
[~,tmp]=NumInt([0;solnew.tau(indU)],[0,prf.sU(end:-1:1)] ,2,true);
solnew.CI_U=tmp;
[~,tmp]=NumInt([0;solnew.tau(indL)],[0,prf.sL]           ,2,true);
solnew.CI_L=tmp;


end


