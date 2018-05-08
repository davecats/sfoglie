function [Tab ] = BlowingComparison(prf,wake,sol,prfB,solB,section,Quantity)
%BLOWINGCOMPARISON gives Plots of Comparison between blowing case and case without blowing
%                  section=0: only gives overview plot (default)
%                  section=1: gives comparison plots for suction side 
%                  section=2: gives comparison plots for pressure side 
%                  Quantity: quantity to be plotted  
%                       'tau'   - All friction related quantities (Cfint,tau, r)  (default)
%                       'Cp'    - Pressure coefficient
%                       'delta' - Boundary layer thickness
%                       'q'     - source terms




if nargin==5; section=0; end
if nargin==6; Quantity='Cf'; end


% reduction in global values


Cd_red=  (solB.Cdrag-sol.Cdrag)/sol.Cdrag ;
Cnu_red=  (solB.Cnu-sol.Cnu)/sol.Cnu ;
CL_red=  (solB.CL-sol.CL)/sol.CL ;

cwp=sol.Cdrag-sol.Cnu;
cwpB=solB.Cdrag-solB.Cnu;
cwp_red=(cwpB-cwp) /cwp;

ratio_red= (solB.CL/solB.Cdrag - sol.CL/sol.Cdrag)*sol.Cdrag/sol.CL;

str={['\Delta C_L=',num2str(100*CL_red),'%'],['\Delta C_d=',num2str(100*Cd_red),'%'],['\Delta C_\nu=',num2str(100*Cnu_red),'%'],['\Delta C_L/C_d=',num2str(100*ratio_red),'%']};


% plot Profile with blowing stuff
xUtr= [prf.nodes.X(sol.iTran(1));prf.nodes.Y(sol.iTran(1))] + sol.tran.Lturb(1)*prf.panels.e(:,sol.iTran(1));
xLtr= [prf.nodes.X(sol.iTran(2));prf.nodes.Y(sol.iTran(2))] - sol.tran.Lturb(2)*prf.panels.e(:,sol.iTran(2)-1);
lls=0.03*prf.c;
endU= xUtr - lls*prf.panels.n(:,sol.iTran(1));
endL= xLtr - lls*prf.panels.n(:,sol.iTran(2)-1);

fig1=PlotProfile( prfB,wake,solB, 4) ;
line([xUtr(1) endU(1)]  , [xUtr(2) endU(2)],'color','k','Linewidth',0.7);
line([xLtr(1) endL(1)]  , [xLtr(2) endL(2)],'color','k','Linewidth',0.7);

text(0.4*prf.c, 0.28*prf.c, str);

% overview table
colHead={'c_a','c_w','c_wR','c_wp','ratio'};
rowHead={'No Blow','Blow','rel. change [%]'};


 TT= [sol.CL , sol.Cdrag , sol.Cnu, cwp ,sol.CL/sol.Cdrag;...
      solB.CL, solB.Cdrag, solB.Cnu,cwpB,solB.CL/solB.Cdrag  ;...
      round(CL_red*100,3), round(Cd_red*100,3),round(Cnu_red*100,3),round(cwp_red*100,3),round(ratio_red*100,3)];

  
Tab=table(TT(:,1),TT(:,2),TT(:,3),TT(:,4),TT(:,5),'RowNames',rowHead,'VariableNames',colHead)

% Tab=uitable;
% Tab.Data=TT;
% Tab.ColumnName=colHead;
% Tab.RowName=rowHead;

if section==0; return; end


% calculate drag reduction coefficients
%--------------------------------------------------------------------

%dS=prfB.sLE-prf.sLE;

%suction side
sU=[0,prf.sU(end:-1:1)];
sBU=[0,prfB.sU(end:-1:1)];

% Interpolate basic and blowing solution to same arclength
B_CFI_U=spline(sBU',[0;solB.tau(prfB.Nle-1:-1:1)],sU');
rU=1-[1; B_CFI_U(2:end)./sol.tau(prf.Nle-1:-1:1)];

B_CFintI_U=spline(sBU',solB.CI_U,sU');
RU=1-B_CFintI_U./sol.CI_U;
RU(1)=0;

%pressure side
sL=[0,prf.sL];
sBL=[0,prfB.sL];

% Interpolate basic and blowing solution to same arclength
B_CFI_L=spline(sBL',[0;solB.tau(prfB.Nle:prf.N)],sL');
rL=1-[1; B_CFI_L(2:end)./sol.tau(prf.Nle:prf.N)];

B_CFintI_L=spline(sBL',solB.CI_L,sL');
RL=1-B_CFintI_L./sol.CI_L;
RL(1)=0;

%
% comparison plots with blowing and without blowing
%---------------------------------------------------

VbU=solB.Vb(1:prfB.Nle-1);
VbL=solB.Vb(prfB.Nle-1:prfB.N);
indB1= find( abs(VbU)>1e-7 );
indB2= find( abs(VbL)>1e-7 );

if strcmp(Quantity ,'Cf') || strcmp(Quantity ,'CF')|| strcmp(Quantity ,'tau')|| strcmp(Quantity ,'Tau')
    
    %-----------  Drag reduction factor ------------------
    figure 
    hold on
    if section==1 % suction side
        plot(sU, rU ,'g')
        plot(sU, RU ,'b')
        if ~isempty(indB1)
            line([prf.sU(indB1(1))   prf.sU(indB1(1))]  , [0  max(rU)+0.03],'color','black','LineStyle','--');
            line([prf.sU(indB1(end)) prf.sU(indB1(end))], [0  max(rU)+0.03],'color','black','LineStyle','--');
        end
        add=' on suction side';
        ylim([min(rU)  max(rU)+0.03])
    elseif section==2 % pressure side
        plot(sL, rL ,'g')
        plot(sL, RL ,'b')
        if ~isempty(indB2)
            line([prf.sL(indB2(1))   prf.sL(indB2(1))]  , [0  max(rL)+0.03],'color','black','LineStyle','--');
            line([prf.sL(indB2(end)) prf.sL(indB2(end))], [0  max(rL)+0.03],'color','black','LineStyle','--');
        end
        add=' on pressure side';
        ylim([min(rL)  max(rL)+0.03])
    end
    title(['Drag-reduction coefficient',add])
    ylabel(' r ') 
    xlabel(' s ')
    legend('local r','global R','blowing region','location','northeast'); 

    %-----------  Integral Cf ------------------
    figure 
    hold on
    if section==1 % suction side
        plot(sU, sol.CI_U ,'g')
        plot(sU, B_CFintI_U,'b')
        if ~isempty(indB1)
            line([prf.sU(indB1(1))   prf.sU(indB1(1))]  , [0  max(sol.CI_U)*1.05],'color','black','LineStyle','--');
            line([prf.sU(indB1(end)) prf.sU(indB1(end))], [0  max(sol.CI_U)*1.05],'color','black','LineStyle','--');
        end
        add=' on suction side';
    elseif section==2 % pressure side
        plot(sL, sol.CI_L ,'g')
        plot(sL, B_CFintI_L,'b')
        if ~isempty(indB2)
            line([prf.sL(indB2(1))   prf.sL(indB2(1))]  , [0  max(sol.CI_L)*1.05],'color','black','LineStyle','--');
            line([prf.sL(indB2(end)) prf.sL(indB2(end))], [0  max(sol.CI_L)*1.05],'color','black','LineStyle','--');
        end
        add=' on pressure side';
    end
    xlim([0,1])
    title(['integral friction coefficient',add])
    ylabel(' [c_f] ') 
    xlabel(' s ')
    legend('without blowing','with blowing','blowing region','location','southeast'); 

    % ------------ tau ----------------

    figure 
    hold on
    plot(prf.nodes.X(1:prf.N), sol.tau(1:prf.N) ,'g')
    plot(prf.nodes.X(1:prf.N), solB.tau(1:prf.N) ,'b')
    title('wall shear stress')
    ylabel(' \tau_w ') 
    xlabel(' x ')
    xlim([0,1])
    legend('without blowing','with blowing','location','northeast'); 

elseif strcmp(Quantity ,'Cp') || strcmp(Quantity ,'CP')
    %               Cp
    %--------------------------------------------------------------------
    figure 
    hold on
    plot(prf.nodes.X(1:prf.N), sol.Cp(1:prf.N) ,'g')
    plot(prf.nodes.X(1:prf.N), solB.Cp(1:prf.N) ,'b')
    title('Pressure coefficient')
    ylabel(' C_p ') 
    xlabel(' x ')
    legend('without blowing','with blowing','location','northeast'); 
    
 elseif strcmp(Quantity ,'delta') 
     % displ and momentum Thickness
     %---------------------------------
     if section==1
        figure 
        hold on
        plot(prf.nodes.X(1:prf.Nle-1), sol.D(1:prf.Nle-1) ,'k')
        plot(prf.nodes.X(1:prf.Nle-1), sol.T(1:prf.Nle-1) ,'b')
        plot(prf.nodes.X(1:prf.Nle-1), solB.D(1:prf.Nle-1) ,'r')
        plot(prf.nodes.X(1:prf.Nle-1), solB.T(1:prf.Nle-1) ,'g')
        if ~isempty(indB1)
            line([prf.xU(indB1(1))   prf.xU(indB1(1))]  , [0  max(solB.D(1:prf.Nle-1))*1.05],'color','black','LineStyle','--');
            line([prf.xU(indB1(end)) prf.xU(indB1(end))], [0  max(solB.D(1:prf.Nle-1))*1.05],'color','black','LineStyle','--');
        end
        title('BL thickness on suction side')
        ylabel(' \delta ') 
        xlabel(' x ')
        legend('$\delta_1$ without blowing','$\delta_2$ without blowing','$\delta_1$ with blowing','$\delta_2$ with blowing','location','northeast'); 
     else
        figure 
        hold on
        plot(prf.nodes.X(prf.Nle-1:prf.N), sol.D(prf.Nle-1:prf.N) ,'k')
        plot(prf.nodes.X(prf.Nle-1:prf.N), sol.T(prf.Nle-1:prf.N) ,'b')
        plot(prf.nodes.X(prf.Nle-1:prf.N), solB.D(prf.Nle-1:prf.N) ,'r')
        plot(prf.nodes.X(prf.Nle-1:prf.N), solB.T(prf.Nle-1:prf.N) ,'g')
        if ~isempty(indB2)
            line([prf.xL(indB2(1))   prf.xL(indB2(1))]  , [0  max(solB.D(prf.Nle:prf.N))*1.05],'color','black','LineStyle','--');
            line([prf.xL(indB2(end)) prf.xL(indB2(end))], [0  max(solB.D(prf.Nle:prf.N))*1.05],'color','black','LineStyle','--');
        end
        title('BL thickness on pressure side')
        ylabel(' \delta ') 
        xlabel(' x ')
        legend('$\delta_1$ without blowing','$\delta_2$ without blowing','$\delta_1$ with blowing','$\delta_2$ with blowing','location','northeast');  
         
         
     end
     
elseif strcmp(Quantity ,'q') || strcmp(Quantity ,'Q') 
    
    % sources
%---------------------------------
delta_q=solB.m-sol.m;
 if section==1
    figure 
    hold on
    plot(prf.nodes.X(1:prf.Nle-1), sol.m(1:prf.Nle-1) ,'k')
    plot(prf.nodes.X(1:prf.Nle-1), solB.m(1:prf.Nle-1),'b')
    plot(prf.nodes.X(1:prf.Nle-1), delta_q(1:prf.Nle-1),'r ')
    if ~isempty(indB1)
        line([prf.xU(indB1(1))   prf.xU(indB1(1))]  , [0  max(solB.m)*1.05],'color','black','LineStyle','--');
        line([prf.xU(indB1(end)) prf.xU(indB1(end))], [0  max(solB.m)*1.05],'color','black','LineStyle','--');
    end
    title('source distribution suction side')
    ylabel(' q ') 
    xlabel(' x ')
    legend('q without blowing','q with blowing','$q_b-q$','location','northwest'); 
    ylim([0  max(solB.m(1:prfB.Nle-1))*1.05])
 else
    figure 
    hold on
    plot(prf.nodes.X(prf.Nle:prf.N), sol.m(prf.Nle:prf.N) ,'k')
    plot(prf.nodes.X(prf.Nle:prf.N), solB.m(prf.Nle:prf.N),'b')
    plot(prf.nodes.X(prf.Nle:prf.N), delta_q(prf.Nle:prf.N),'r ')
    if ~isempty(indB2)
        line([prf.xU(indB2(1))   prf.xU(indB2(1))]  , [0  max(solB.m)*1.05],'color','black','LineStyle','--');
        line([prf.xU(indB2(end)) prf.xU(indB2(end))], [0  max(solB.m)*1.05],'color','black','LineStyle','--');
    end
    title('source distribution pressure side')
    ylabel(' q ') 
    xlabel(' x ')
    legend('q without blowing','q with blowing','$q_b-q$','location','northwest'); 
    ylim([0  max(solB.m(prfB.Nle:prf.N))*1.05]) 
     
 end
    
% 
% figure 
% hold on
% plot(sges, sol.m ,'k')
% plot(sges, solB.m,'b')
% plot(sges, delta_m,'r --')
% ylim([0  max(solB.q)*1.05])
% legend('q without blowing','q with blowing','$q_b-q$','location','northwest')


end

