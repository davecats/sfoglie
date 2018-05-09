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
%                       'H12'   - shape parameter
%                       'H32'   - shape parameter
%                       'CD'    - dissipation coefficient
%                       'D'     - dissipation integral
%                       'ReT'   - Reynoldsnumber based on momentum thickness



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


legendstr={'without blowing','with blowing'};

if strcmp(Quantity ,'Cf') || strcmp(Quantity ,'CF')|| strcmp(Quantity ,'tau')|| strcmp(Quantity ,'Tau')
    titlestr='wall shear stress';
    ystr=' C_f= \tau_w / \rho U_\infty^2 ';
    yNo=sol.tau;
    yBlow=solB.tau;
elseif strcmp(Quantity ,'Cp') || strcmp(Quantity ,'CP')
    titlestr='pressure coefficient';
    ystr=' C_p ';
    yNo=sol.Cp;
    yBlow=solB.Cp;
elseif strcmp(Quantity ,'H12') 
    titlestr='Shape parameter';
    ystr=' H_1_2 ';
    yNo=sol.HK;
    yBlow=solB.HK;
elseif strcmp(Quantity ,'H12') 
    titlestr='Shape parameter';
    ystr=' H_3_2 ';
    yNo=sol.HS;
    yBlow=solB.HS;
elseif strcmp(Quantity ,'CD') || strcmp(Quantity ,'cD')
    titlestr='Dissipation coefficient';
    ystr=' C_D ';
    yNo=sol.CD;
    yBlow=solB.CD;
elseif strcmp(Quantity ,'D')
    titlestr='Dissipation integral';
    ystr=' D ';
    yNo=sol.CD.*sol.U.^3;
    yBlow=solB.CD.*solB.U.^3;
elseif strcmp(Quantity ,'delta')
    titlestr='Boundary Layer thickness';
    ystr=' \delta ';
    yNo=sol.D;
    yBlow=solB.D;
    
    y2No=sol.T;
    y2Blow=solB.T;
    legendstr={'$\delta_1$ without blowing','$\delta_1$ with blowing','$\delta_2$ without blowing','$\delta_2$ with blowing'};
elseif strcmp(Quantity ,'q') || strcmp(Quantity ,'Q')
    titlestr='source contribution';
    ystr=' q ';
    yNo=sol.m;
    yBlow=solB.m;  
elseif strcmp(Quantity ,'Ret') || strcmp(Quantity ,'ReT')
    titlestr='Reynoldsnumber based on momentum thickness';
    ystr=' Re_t';
    yNo=sol.Ret;
    yBlow=solB.Ret;    
end

if section==1
    ind=1:prf.Nle-1;
elseif section==2
    ind=prf.Nle:prf.N;
end
xstr='x';

%upper limit
upl= max([yNo(ind);yBlow(ind)])*1.05;
%lower limit
lowl= min( 0,min([yNo(ind);yBlow(ind)]) );




% plot
figure 
hold on
title(titlestr) 
xlabel(xstr)
ylabel(ystr)
plot(prf.nodes.X(ind),yNo(ind), 'b')
plot(prf.nodes.X(ind),yBlow(ind), 'b --')
if strcmp(Quantity ,'delta')
    plot(prf.nodes.X(ind),y2No(ind),'r')
    plot(prf.nodes.X(ind),y2Blow(ind),'r --')  
elseif strcmp(Quantity ,'q') || strcmp(Quantity ,'Q')
    delta_q=solB.m-sol.m;
    legendstr={'q without blowing','q with blowing','$q_b-q$'};
    plot(prf.nodes.X(ind),delta_q(ind), 'r')
end
% blowing region
if ~isempty(indB1) && section==1
    line([prf.xU(indB1(1))   prf.xU(indB1(1))]  , [lowl  upl],'color','black','LineStyle','--');
    line([prf.xU(indB1(end)) prf.xU(indB1(end))], [lowl  upl],'color','black','LineStyle','--');
end
if ~isempty(indB2) && section==2
    line([prf.xL(indB2(1))   prf.xL(indB2(1))]  , [lowl  upl],'color','black','LineStyle','--');
    line([prf.xL(indB2(end)) prf.xL(indB2(end))], [lowl  upl],'color','black','LineStyle','--');
end
%legendstr={legendstr,'blowing region'};
legend(legendstr,'location','best'); 
xlim([0,max(prf.nodes.X)])


if strcmp(Quantity ,'Cp') || strcmp(Quantity ,'CP')
    ind=1:prf.N;
    figure 
    hold on
    title(titlestr) 
    xlabel(xstr)
    ylabel(ystr)
    plot(prf.nodes.X(ind),yNo(ind), 'b')
    plot(prf.nodes.X(ind),yBlow(ind), 'b --')
    legend(legendstr,'location','best'); 
    xlim([0,max(prf.nodes.X)])  
end



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

end

 
end
    
% 
% figure 
% hold on
% plot(sges, sol.m ,'k')
% plot(sges, solB.m,'b')
% plot(sges, delta_m,'r --')
% ylim([0  max(solB.q)*1.05])
% legend('q without blowing','q with blowing','$q_b-q$','location','northwest')



