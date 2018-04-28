init
parameters
% inviscid solution 
[InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);

% InvRef is a struct containing all information from inviscid solution 
    %-> can be committed by next use of airfoil function as input to get shorter calculation process
[sol,prf,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);


%---------------------------------------------------------

% set reference case
CL_ref=sol.CL; CD_ref=sol.Cdrag; tr_ref=sol.tran.x;

% use tripping at transition location ofuncontrolled case (not necessary)
solRef=sol;prfRef=prf;
tri.active=[true;true];
tri.x=solRef.tran.x;
flo.nkrit=9;

% blo.x={0.048402*prf.c;         
%        0.050529*prf.c;};
% blo.active=true;
% [solB,prfB,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);
% BlowingComparison(prf,flo.wake,sol,prfB,solB,1);

%
%% particel swarm optimization
%--------------------------------------------------------------
% inviscid solution -> only calculate once because it does not change
[InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);


% fixed parameters for optimization
global OptiInputStruct % musst be global to be avaivable in Optimization function handle
OptiInputStruct=InvRef;
OptiInputStruct.blo=blo;OptiInputStruct.tri=tri;OptiInputStruct.eng=eng;

global PSOerg
PSOerg=[];
handle=@fOPT;
lowerConstraint=-0.1*[1,1];
upperConstraint=1.1*[1,1];
PS_Options=optimoptions('particleswarm','SwarmSize',20,'PlotFcn','pswplotbestf','OutputFcn',@PSOout);

% optimization process
rng default 
[p,val,exitCond,out]=particleswarm(handle,2,lowerConstraint,upperConstraint,PS_Options);



%   Plot
%
%scatter3(PSOerg(:,1),PSOerg(:,2),PSOerg(:,3))

% Plot of optimization parameters
[xm,ym]=meshgrid(-0.1:0.02:1.1);
zm=griddata(PSOerg(:,1),PSOerg(:,2),PSOerg(:,3),xm,ym);
figure
hold on
contourf(xm,ym,zm)
scatter(PSOerg(:,1),PSOerg(:,2),'k x')

% compare optimal case
blo.x={p(1)*prf.c;         
       p(2)*prf.c;};
blo.active=true;
[sol,prf,flo,~,~,~]=airfoil(prf,flo,tri,blo,eng);
BlowingComparison(prfRef,flo.wake,solRef,prf,sol,1);


%dlmwrite('./optimization/Re4e5_alfa12_Tr_NoSkw.txt',PSOerg)

%polarF = './optimization/Re4e5_Polar_Tr_NoSkw.txt';
%dlmwrite(polarF,[12, p,sol.CL,sol.Cdrag,sol.Cnu,sol.CL/sol.Cdrag], '-append')


% %%
% % % Plots
% %--------------------------------------------------
% 
% % quantities that can be plotted: 
% %'tau','Cf','Cfint','Cp','delta' ,'U', 'CD','D','H12' ,'H32', 'Ret'
% % section: 1 - pressure side, 2 - suction side, 3 - wake, 4 - suction and pressure side , 5 - all sections
% % OverArclength=true ->  plot over arclength s instead of x                 
% % section=4;
% % OverArclength=false;
% % PlotStuff(prf,flo.wake,sol, 'tau',section,OverArclength);
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
% % mode 1: only plots overview with reduction in CL and Cd
% % mode 2: plots overview + plots for suction side
% % mode 3: plots overview + plots for pressure side
% mode=1;
% BlowingComparison(prf,flo.wake,sol,prfB,solB,1);
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














