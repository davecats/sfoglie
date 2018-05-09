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
% tri.active=[true;true];
% tri.x=solRef.tran.x;
% flo.nkrit=9;



%% particel swarm optimization
%--------------------------------------------------------------
% 
%% one blowing region on suction and pressure side -> position optimized

% % inviscid solution -> only calculate once because it does not change
% [InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);
% 
% 
% % fixed parameters for optimization
% global OptiInputStruct % musst be global to be avaivable in Optimization function handle
% OptiInputStruct=InvRef;
% OptiInputStruct.blo=blo;OptiInputStruct.tri=tri;OptiInputStruct.eng=eng;
% 
% % variing the position of two blowing regions on suction and pressure side
% global PSOerg
% PSOerg=[];
% handle=@fOPT;
% lowerConstraint=-0.1*[1,1];
% upperConstraint=1.1*[1,1];
% particlenumber=20;    % number of particle that are used
% PS_Options=optimoptions('particleswarm','SwarmSize',particlenumber,'PlotFcn','pswplotbestf','OutputFcn',@PSOout);
% 
% % optimization process
% rng default 
% [p,val,exitCond,out]=particleswarm(handle,2,lowerConstraint,upperConstraint,PS_Options);
% 
% 
% 
% %   Plot
% % 
% %scatter3(PSOerg(:,1),PSOerg(:,2),PSOerg(:,3))
% ind=find(PSOerg(:,3)==0);
% PSOerg(ind,:)=[];
% % Plot of optimization parameters
% [xm,ym]=meshgrid(-0.1:0.02:1.1);
% zm=griddata(PSOerg(:,1),PSOerg(:,2),PSOerg(:,3),xm,ym);
% figure
% hold on
% contourf(xm,ym,zm)
% scatter(PSOerg(:,1),PSOerg(:,2),'k x')
% 
% % compare optimal case
% blo.x={p(1)*prf.c;         
%        p(2)*prf.c;};
% blo.active=true;
% [sol,prf,flo,~,~,~]=airfoil(prf,flo,tri,blo,eng);
% BlowingComparison(prfRef,flo.wake,solRef,prf,sol,2,'delta');
% 
% % write out values
% dlmwrite('./optimization/Re4e5_alfa16_N15_NoSkw.txt',PSOerg)
% 
% polarF = './optimization/Re4e5_Polar_N15_NoSkw.txt';
% dlmwrite(polarF,[16, p,sol.CL,sol.Cdrag,sol.Cnu,sol.CL/sol.Cdrag], '-append')




%% one blowing region only on suction side-> position, length, intensity optimized
% 
% % inviscid solution -> only calculate once because it does not change
% [InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);
% % fixed parameters for optimization
% global OptiInputStruct % musst be global to be avaivable in Optimization function handle
% OptiInputStruct=InvRef;
% OptiInputStruct.blo=blo;OptiInputStruct.tri=tri;OptiInputStruct.eng=eng;
% 
% % % variing the position, length and intensity of a blowing region on the suction side
% global PSOerg
% PSOerg=[];
% handle=@fOPT2;
% lowerConstraint=[0,0,-0.01];
% upperConstraint=[1,1,0.01];
% particlenumber=80;%300;    % number of particle that are used
% PS_Options=optimoptions('particleswarm','SwarmSize',particlenumber,'PlotFcn','pswplotbestf','OutputFcn',@PSOout);
% 
% % optimization process
% rng default 
% [p,val,exitCond,out]=particleswarm(handle,3,lowerConstraint,upperConstraint,PS_Options);
% 
% % delete dummy values in regions out of the constraint x_s + L_b <1
% ind=find(PSOerg(:,4)==0);
% PSOerg(ind,:)=[];
% 
% % % write out
% % dlmwrite('./optimization/Re4e5_alfa2_Tr_NoSkw_suction.txt',PSOerg)
% % 
% % polarF = './optimization/Re4e5_Polar_Tr_NoSkw_suction.txt';
% % dlmwrite(polarF,[2, p,sol.CL,sol.Cdrag,sol.Cnu,sol.CL/sol.Cdrag], '-append')
% 
% % compare optimum case
% xm=p(1)+p(2)/2;
% blo.active=true;                 
% blo.L= {[p(2)]*OptiInputStruct.prf.c;              %  length of blowing area: suction side
%         [0.1]*OptiInputStruct.prf.c;};            %  length of blowing area: pressure side
% blo.x= {xm*OptiInputStruct.prf.c;          %  midpoint of blowing area
%         0.1*OptiInputStruct.prf.c;};             
% blo.A= {[p(3)]*OptiInputStruct.flo.Uinfty;
%         [0.0]*OptiInputStruct.flo.Uinfty};
% 
% [solB,prfB,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);
% BlowingComparison(prf,flo.wake,sol,prfB,solB,1);
% 
% 
% %% Plots
% 
% figure
% scatter3(PSOerg(:,1),PSOerg(:,2),PSOerg(:,3),'.')
% % Plot of optimization parameters
% pp1=0:0.02:1;
% pp2=0:0.02:1;in=find(pp2>p(2),1,'first'); pp2=[pp2(1:in-1),p(2),pp2(in:end)];
% pp3=-0.01:0.0004:0.01;
% [xm,ym,zm]=meshgrid(pp1,pp2,pp3);
% erg=griddata(PSOerg(:,1),PSOerg(:,2),PSOerg(:,3),PSOerg(:,4),xm,ym,zm);
% 
% % Iso-surfaces where optimum point lays in
% % x_s= konst
% x_iso1=squeeze(ym(:,1,:));
% y_iso1=squeeze(zm(:,1,:));
% z_iso1=squeeze(erg(:,1,:));
% 
% % L_b= konst
% x_iso2=squeeze(xm(5,:,:));
% y_iso2=squeeze(zm(5,:,:));
% z_iso2=squeeze(erg(5,:,:));
% 
% % v_w= konst
% x_iso3=squeeze(xm(:,:,end));
% y_iso3=squeeze(ym(:,:,end));
% z_iso3=squeeze(erg(:,:,end));
% 
% 
% gg= 1-x_iso3 - y_iso3;
% z_iso3(gg<1e-7)=NaN;
% 
% % plot of iso surfaces
% figure
% hold on
% contourf(x_iso1,y_iso1,z_iso1)
% plot(p(2),p(3),'x')
% xlabel('length of blowing region')
% ylabel('intensity v_w')
% 
% figure
% hold on
% contourf(x_iso2,y_iso2,z_iso2)
% plot(p(1),p(3),'x')
% xlabel('startinpoint of blowing region')
% ylabel('intensity v_w')
% 
% figure
% hold on
% contourf(x_iso3,y_iso3,z_iso3)
% plot(p(1),p(2),'x')
% xlabel('startinpoint of blowing region')
% ylabel('length of blowing region')



%% two control regions with arbitrary positon -> position, length, intensity optimized

% inviscid solution -> only calculate once because it does not change
[InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);
% fixed parameters for optimization
global OptiInputStruct % musst be global to be avaivable in Optimization function handle
OptiInputStruct=InvRef;
OptiInputStruct.blo=blo;OptiInputStruct.tri=tri;OptiInputStruct.eng=eng;

% % variing the position, length and intensity of a blowing region on the suction side
global PSOerg
PSOerg=[];
handle=@fOPT3;
lowerConstraint=[0,0,0,0,-0.01,-0.01];
upperConstraint=[1,1,1,1, 0.01, 0.01];
particlenumber=450;    % number of particle that are used
PS_Options=optimoptions('particleswarm','SwarmSize',particlenumber,'PlotFcn','pswplotbestf','OutputFcn',@PSOout);

% optimization process
rng default 
[p,val,exitCond,out]=particleswarm(handle,6,lowerConstraint,upperConstraint,PS_Options);

 ind=find(PSOerg(:,4)==0);
 PSOerg(ind,:)=[];

% write out
dlmwrite('./optimization/Op3_alfa2_Tr_NoSkw_suction.txt',PSOerg)

blo.active=true;
blo.ArcLengthMode=true;
% arclengths in percent of smax
blo.s_change=        [p(1),p(2),p(3),p(4)]; 
% intensitie for s_i< s < s_i+1
blo.NewA= OptiInputStruct.flo.Uinfty*[p(5),0,p(6),0];

[solB,prfB,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);

BlowingComparison(prf,flo.wake,sol,prfB,solB,2,'tau');













