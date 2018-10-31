% 
%      //////  /////  ////   /////   //     //  /////
%     //   // //    //   // //  //  //     //  //
%     //     //    //   // //      //     //  //
%      //   ////  //   // // //// //     //  /////
%  //   // //    //   // //   // //     //  //
%   ///// //      ////    ////  ////// //  /////
%
%   
% This program computes the incompressible flow around an airfoil 
% with and without uniform blowing
%
% ------------------------------------------------------------------------------
% If you use this code and find it helpful please cite
% 
% M.Reder, A. Stroh and D. Gatti, "Preliminary study of flow control via uniform 
% blowing on airfoils with a boundary element method", Notes on Numerical Fluid 
% Mechanics and Multidisciplinary Design, (submitted, 2018)
% ------------------------------------------------------------------------------
% 
% Copyright (C) 2018  Dr. Davide Gatti & M. Reder
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Contacts:
%
% Dr. Davide Gatti  
% davide.gatti [at] kit.edu
% msc.davide.gatti [at] gmail.com
% 
% Karlsruhe Institute of Technology
% Institute of Fluid Dynamics
% Kaiserstraße 10
% 76131 Karlsruhe 
% 


init        % initialize: add paths etc.
parameters  % load parameters for calculation
% inviscid solution 
[InvRef.prf,InvRef.flo,InvRef.CoeffMatrix,InvRef.Uinv,InvRef.sges ] = InviscidSolution(prf,flo,eng);

% InvRef is a struct containing all information from inviscid solution 
    %-> can be committed by next use of airfoil function as input to get shorter calculation process
[sol,prf,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);

%---------------------------------------------------------

% set reference case without uniform blowing
CL_ref=sol.CL; CD_ref=sol.Cdrag; tr_ref=sol.tran.x;ratio_ref=CL_ref/CD_ref;

% use tripping at transition location of uncontrolled case 
%   ->neglect effect of blowing on transition, uncomment to consider transition changes
solRef=sol;prfRef=prf;
tri.active=[true;true];
tri.x=solRef.tran.x;
flo.nkrit=9;

% turn on uniform blowing
blo.active=true;                     
[solB,prfB,flo,~,~,~,~]=airfoil(prf,flo,tri,blo,eng,InvRef);

% compare
BlowingComparison(prf,flo.wake,sol,prfB,solB,1,'delta');

% %% Plotting Functions avaivable 
%
% % % plot quantities
% %------------------------------------------------------------------
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
% % Overview-plots of profiles
% %------------------------------------------------------------------
% % mode 1: plots Profil and wake with displacement thickness and shows CL, Cd usw.
% % mode 2: plots Profil and wake with the node positions used
% % mode 3: plots Profil with transition locations and blowing distribution in blowing case
% %------------------------------------------------------------------
% mode=1;
% PlotProfile(prf,flo.wake,sol, mode);
%
% 
% % Comparison between blowing case and reference case without blowing
% %------------------------------------------------------------------
% % section 0: only plots overview with reduction in CL and Cd
% % section 1: plots overview + plots for suction side
% % section 2: plots overview + plots for pressure side
% % Quantities that can be plotted: 'Cf','Cp','delta','q','H12','H32','CD','D','Ret'  
% %------------------------------------------------------------------
% section=1;Quantity='Cf'; 
% BlowingComparison(prf,flo.wake,sol,prfB,solB,section,Quantity);
% 
% % Plot line forces 
% %------------------------------------------------------------------
% % mode=0: sum of lower and upper airfoil part
% % mode=1: upper airfoil part
% % mode=2: lower airfoil part
% % if sol2 comitted -> plot both solutions
% %------------------------------------------------------------------
% PlotForces(prf,flo,sol,mode, solRef)
%
% % Plot stress/force vectors
% %------------------------------------------------------------------
% % each_nth: only plot each n-th node (default 1 -> every node)
% % total   =true: plot the total stresses consisting of pressure and shear stresses 
% % pressure=true: plot the pressure stresses
% % shear   =true: plot the shear stresses
% % force = true : plots force instead of stress vectors (integrated over area)
% %------------------------------------------------------------------
% each_nth=3,total=true; pressure=false; shear=false; force=false;
% PlotStressVectors( prf,sol,flo,each_nth,total, pressure, shear, force )



return
%% Calculation of polar curves
%---------------------------------------------------------
% loop over alfas + write out  for polar curves
step=1;start=0;
for k=start:step:12

    flo.alfa= k *pi/180;

    disp(['ALFA=',num2str(k)])
    disp('--------------------------------------------')
    
    
    % evaluation
    [sol,prf,flo,~,~,~]=airfoil(prf,flo,tri,blo,eng);

    if k==start
        CP=sol.Cp(1:prf.N);
        TAU=sol.tau(1:prf.N);
        CL   =sol.CL ;
        Cnu  =sol.Cnu;
        Cdrag=sol.Cdrag;
        Cd_p=sol.Cdrag-sol.Cnu;
        alph=k;
        
        xTT_up= [sol.tran.x(1), sol.xseparation(1), sol.xreattach(1) ];
        xTT_lo= [sol.tran.x(2), sol.xseparation(2), sol.xreattach(2) ];
    else
        CP=[CP,sol.Cp(1:prf.N)];
        TAU=[TAU,sol.tau(1:prf.N)];
        CL=[CL;sol.CL];
        Cnu=[Cnu;sol.Cnu];
        Cdrag=[Cdrag;sol.Cdrag];
        Cd_p=[Cd_p;sol.Cdrag-sol.Cnu];
        alph=[alph;k];
        
        xTT_up=[xTT_up;sol.tran.x(1), sol.xseparation(1), sol.xreattach(1) ];
        xTT_lo=[xTT_lo;sol.tran.x(2), sol.xseparation(2), sol.xreattach(2) ];

    end

end

%% Even distributed paremeters for blowing location  
%
% N=20;
% xTop=linspace(-0.1,1,N);
% xBot=linspace(-0.1,1,N);
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
% xlabel('x^{TOP}')
% ylabel('x^{BOTTOM}')














