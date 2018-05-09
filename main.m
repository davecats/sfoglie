init
parameters
% inviscid solution 
[InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);

% InvRef is a struct containing all information from inviscid solution 
    %-> can be committed by next use of airfoil function as input to get shorter calculation process
[sol,prf,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);

%---------------------------------------------------------

% set reference case
CL_ref=sol.CL; CD_ref=sol.Cdrag; tr_ref=sol.tran.x;ratio_ref=CL_ref/CD_ref;

% use tripping at transition location ofuncontrolled case (not necessary)
solRef=sol;prfRef=prf;
tri.active=[true;true];
tri.x=solRef.tran.x;
flo.nkrit=9;

%p=[0,0.08016,0.01];
p=[0.8292,0.1262,-0.01];
xm=p(1)+p(2)/2;
blo.active=true;                 
blo.L= {[p(2)]*prf.c;              %  length of blowing area: suction side
        [0.1]*prf.c;};            %  length of blowing area: pressure side
blo.x= {xm*prf.c;          %  midpoint of blowing area
        0.1*prf.c;};             
blo.A= {[p(3)]*flo.Uinfty;
        [0.0]*flo.Uinfty};
    
[solB,prfB,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);

BlowingComparison(prf,flo.wake,sol,prfB,solB,1,'delta');

% %%
% % % Plots
% %--------------------------------------------------
% 
% % quantities that can be plotted: 
% %'tau','Cf','Cfint','Cp','delta' ,'U', 'CD','D','H12' ,'H32', 'Ret'
% % section: 1 - suction side, 2 - pressure side, 3 - wake, 4 - suction and pressure side , 5 - all sections
% % OverArclength=true ->  plot over arclength s instead of x                 
% % section=4;
% % OverArclength=false;
% % Quantity='tau';  
% % PlotStuff(prf,flo.wake,sol, Quantity,section,OverArclength);
% %------------------------------------------------------------------
% 
% % default: section = 4, OverArclength=false;
% PlotStuff(prf,flo.wake,sol, 'tau');
% 
% 
% % Plots of Profiles
% % mode 1: plots Profil and wake with displacement thickness and shows CL, Cd usw.
% % mode 2: plots Profil and wake witch the nodepositions used
% % mode 3: plots Profil with transition locations andblowing distribution in blowing case
% mode=1;
% PlotProfile(prf,flo.wake,sol, mode);
% 
% % Comparison between blowing case and reference case without blowing
% % section 0: only plots overview with reduction in CL and Cd
% % section 1: plots overview + plots for suction side
% % section 2: plots overview + plots for pressure side
% % Quantities that can be plotted: 'Cf','Cp','delta','q','H12','H32','CD','D','ReT'  
% section=1;Quantity='Cf'; 
% BlowingComparison(prf,flo.wake,sol,prfB,solB,section,Quantity);
% 


return
%% Calculation of polar curves
%---------------------------------------------------------
% loop over alfas + write out  for polar curves
step=1;
for k=0:step:16

    flo.alfa= k *pi/180;

    aTMP=round(flo.alfa*180/pi,1);
    disp(['ALFA=',num2str(aTMP)])
    disp('--------------------------------------------')
    
    
    % evaluation
    [sol,prf,flo,~,~,~]=airfoil(prf,flo,tri,blo,eng);

    if k==0
        CP=sol.Cp(1:prf.N);
        TAU=sol.tau(1:prf.N);
        CL   =sol.CL ;
        Cnu  =sol.Cnu;
        Cdrag=sol.Cdrag;
        alph=aTMP;
        
        xTT= [sol.tran.x(1), sol.xseparation(1), sol.xreattach(1) ];
        
        TST= DragCoeff(sol.T(end),sol.HK(end),sol.U(end),1);
    else
        CP=[CP,sol.Cp(1:prf.N)];
        TAU=[TAU,sol.tau(1:prf.N)];
        CL=[CL;sol.CL];
        Cnu=[Cnu;sol.Cnu];
        Cdrag=[Cdrag;sol.Cdrag];
        alph=[alph;aTMP];
        
        xTT=[xTT;sol.tran.x(1), sol.xseparation(1), sol.xreattach(1) ];
        
        TST=[TST; DragCoeff(sol.T(end),sol.HK(end),sol.U(end),1)];
    end

end

%% Even distributed paremeters for blowing location  
%
% N=20;
% xTop=linspace(-0.05,1.05,N);
% xBot=linspace(-0.05,1.05,N);
% 
% tri.x=tr_ref;
% tri.active=[true; true];
% blo.active=true;
% 
% for i=1:N
%     for j=1:N
%         blo.x= {[xTop(i)]*prf.c;          %  midpoint of blowing area
%                 [xBot(j)]*prf.c;};  
%         [sol,~,~,~,~,~]=airfoil(prf,flo,tri,blo,eng);
%         CL(i,j)=sol.CL;
%         CD(i,j)=sol.Cdrag;
%     end
% end
% 
% contourf(xTop,xBot,(real(CL(:,:))./CD)'/(CL_ref/CD_ref))
% xlabel('x_M^{TOP}')
% ylabel('x_M^{BOTTOM}')














