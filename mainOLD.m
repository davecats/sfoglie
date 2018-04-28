 % Script for airfoil calculation with plots and a comparison blowing/ no blowing 

%init
%%
%  air foil geometry and flow Parameters
%------------------------------------------------------------------------------
% BlowingComparison(prfE,flo.wake,sol,prfB,solB,1);
% 

%  prf and panels
%  ------------------
prf.naca = [4 4 1 2];        %  NACA 4-digit prf 
prf.noSkew  = true;          %  if true neglects prf skewness
prf.sharpTE = false;         %  if true modifies NACA prf for shart trailing edge
prf.c = 1;                   %  prf chord length 
prf.M = 180;                 %  number of control points for each surface
prf.pmode = 1;               %  node distribution mode: 1 (more nodes in middle) 2 (more nodes at LE and TE)

%  Flow
%  ----
flo.alfa = 2*pi/180;         %  Angle of attack 
flo.invisc = false;          %  only inviscid solution 
flo.Uinfty=1;                %  velocity of outer flow
flo.Re= 4e5;                 %  Chord Reynoldsnumber
flo.nkrit=  0.15;%          %  critical amplification exponent for transition

%  Tripping
%  --------
tri.active=[ false;...         %  tripping on suction side 
             false];           %  tripping on pressure side 

tri.x = [ 0.145;...            %  tripping location on suctoin side
          0.29 ]*prf.c;        %  tripping location on pressure side


% ----------- blowing -----------------------
withBlowing=[true;...  % blowing on suction side  
             false];   % blowing on pressure side 
% blowing region
% startpoint    
xBstart= [0.25;...
          0.9]* prf.c;
% end point     
xBend  = [0.5;...
          1]* prf.c;
      
% blowing intensity      
intensity=[0.001;...
           0.001]* flo.Uinfty;

pressureCor=false;%true; %  include correction Term for pressure

%  Newton Solver
%  --------------
eng.it=25;                   % maximum number of Newton step iterations
eng.tranEQ=2;                % 1) xfoil transition EQ 2nd Order 2) xfoil transition EQ 1nd Order 3) modified transition EQ
eng.tol=5e-4;                % tolerance of Newton method


%----------------------------------------------------------------------------
%  inviscid solution
%----------------------------------------------------------------------------

[prf,flo,CoeffMatrix,Uinv,sges ] = InviscidSolution(prf,flo,eng);

% Plot inviscid pressure coefficient
if flo.invisc
    Cp=1-flo.gamma.^2;
    CL=getCL(prf,flo.gamma,flo.alfa,flo.Uinfty);
    if prf.sharpTE; Cp=[Cp; Cp(1)]; end
    figure()
    hold on; box on
    plot([prf.xU,prf.xL(1)], Cp(1:Nle));
    plot(prf.xL,  Cp(Nle:end));
    xlim([0 1]);
    legend('$C_p=1-\gamma(s)^2$ Saugseite','$C_p=1-\gamma(s)^2$ Druckseite');
    
   figure; 
   hold on; box on;
   plot([prf.panels.X]',[prf.panels.Y]','k','Linewidth',2);
   plot(flo.wake.x,flo.wake.y,'b');
   axis equal; xlabel('x'); ylabel('y') 
   return;
end

% [L,U,ind] = lu(flo.Ages,'vector');
% ALU= U + L-eye(size(flo.Ages));


%----------------------------------------------------------------------------
%  viscous solution
%----------------------------------------------------------------------------


%   without blowing
%----------------------------------------------

% initial solution of Boundary Layer for global Newton Method
Vb=zeros(size(Uinv));

% initial solution guess
ini = GetInitialSolution( prf,flo, tri, eng, Uinv, Vb, 2);
%  coupled boundary layer and potential flow solution
[sol, prfE]=NewtonEq( prf,flo,eng,ini,CoeffMatrix.D,Uinv,pressureCor);

% If no convergence -> try with different approach for Transition panel EQ
if sol.residual>eng.tol
    eng.tranEQ=true;
    disp('Not converged, Try different Transition approach ')
    disp('------------------------------------------------ ')
    if sol.residual<1
        [solTST, prfTST]=NewtonEq( prfE,flo,eng,sol,CoeffMatrix.D,Uinv,pressureCor);
    else
        [solTST, prfTST]=NewtonEq( prf,flo,eng,ini,CoeffMatrix.D,Uinv,pressureCor);
    end
    
    if solTST.residual < sol.residual
       sol=solTST; prfE=prfTST;
    end
    clear solTST  prfTST TranEQ2
end
%--------------------------------------------------------------


% plot

% PlotStuff(prfE,flo.wake,sol, 'tau');
% PlotStuff(prfE,flo.wake,sol, 'Cp');

% PlotProfile(prfE,flo.wake,sol, 2);


%%

%   with blowing
%----------------------------------------------


%  tri.active=[ true; true];  
%  tri.x=sol.tran.x;
%  flo.nkrit=  9;% 

if ~withBlowing(1) && ~withBlowing(2); return; end

disp('With Blowing')
disp('-----------------------')

Vb=zeros(size(Uinv));

% get nodes with blowing and set the blowing velocity vector
if withBlowing(1)
    indB1= find( prf.xU < xBend(1) & prf.xU > xBstart(1));
    Vb(indB1)=intensity(1);
end
if withBlowing(2)
    indB2= find( prf.xL < xBend(2) & prf.xL > xBstart(2));
    indB2= indB2 + (prf.Nle-1)*ones(size(indB2));
    Vb(indB2)=intensity(2);
end

% initial solution guess
iniB = GetInitialSolution( prf,flo, tri, eng, Uinv, Vb, 2);
%  coupled boundary layer and potential flow solution
[solB, prfB]=NewtonEq( prf,flo,eng,iniB,CoeffMatrix.D,Uinv,pressureCor);

%--------------------------------------------------------------

% PlotStuff(prfB,flo.wake,solB, 'tau');
% PlotStuff(prfB,flo.wake,solB, 'Cp');

% PlotProfile(prfB,flo.wake,solB, 3);



% % Plots with comparison blowing/ no blowing
% BlowingComparison(prfE,flo.wake,sol,prfB,solB,1);


% % calculate drag reduction coefficients
% %--------------------------------------------------------------------
% 
% %suction side
% sU=[0,prfE.sU(end:-1:1)];
% sBU=[0,prfB.sU(end:-1:1)];
% 
% % Interpolate basic and blowing solution to same arclength
% B_CFI=spline(sBU',[0;solB.tau(prfB.Nle-1:-1:1)],sU');
% rU=1-[1; B_CFI(2:end)./sol.tau(prfE.Nle-1:-1:1)];
% 
% B_CFintI=spline(sBU',solB.CI_U,sU');
% RU=1-B_CFintI./sol.CI_U;
% RU(1)=0;


%%

% Write out 

writeOut=false;

if writeOut
    
    tmp=flo.Re; i=0;
    while tmp>1
        tmp=tmp/10;
        i=i+1;  
    end
    first=round(tmp*10);
    str1=['Re',num2str(first),'e',num2str(i-1)];
    alf=round(flo.alfa*180/pi);
    str2=['_alfa',num2str(alf)];

    if trip(1)
        tr=round(xtrip(1)*100);
        stmp=num2str(tr);
        if strcmp(stmp(2),'0'); stmp=stmp(1); end
        str3=['_Trip0',stmp];
    else
        if nkrit>1
            str3=['_N',num2str(round(nkrit))];
        else
            tmp=100*flo.nkrit;
            stmp=num2str(tmp);
            if strcmp(stmp(2),'0'); stmp=stmp(1); end
            str3=['_N0',stmp];
        end

    end

    str=[str1,str2,str3];

    Xges=[prf.nodes.X'; flo.wake.x];
    Yges=[prf.nodes.Y'; flo.wake.y];

    datBlow  = [Xges,Yges,sges', solB.D,solB.T,solB.U,solB.Cp,solB.tau, solB.Cf ];
    datNoBlow= [Xges,Yges,sges', sol.D,sol.T,sol.U,sol.Cp,sol.tau, sol.Cf ];

     dlmwrite(['./ComparisonData/Blowing_',str], datBlow);
     dlmwrite(['./ComparisonData/NoBlowing_',str], datNoBlow);
end
%--------------------------------------------------------------------






