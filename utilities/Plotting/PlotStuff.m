function [  ] = PlotStuff( prf,wake,sol ,Quantity , section, OverArclength )
%PLOTSTUFF  Function for Plotting 
%           section: 1 - pressure side
%                    2 - suction side
%                    3 - wake
%                    4 - suction and pressure side (default)
%                    5 - all sections
%           OverArclength: if true plots over arclength instead of x
%           Quantity: quantity to be plotted    
%                'tau'   - wall shearstress normalized with rho*Uinfty^2
%                'Cf'    - wall friction coefficient normalized with 0.5*rho*U^2  
%                'Cfint' - integral wall friction coefficient int_0^x cf dx
%                'Cp'    - pressurecoefficient
%                'delta' - diplacement and momentum thickness
%                'U'     - boundary layer edge velocity   
%                'CD'    - Dissipation coefficient
%                'D'     - Dissipation integral   
%                'H12'   - shapeparameter H12
%                'H32'   - shapeparameter H32
%                'Ret'   - Reynoldsnumber based on U and momentumthickness
%                'q'     - sourceterms
%                'm'     - mass defect   
%                'n'     - amplification exponent
%                'ctau'  - shear stress coefficient

%   total arclength vector
if nargin<6
    OverArclength=false;
end
if nargin<5
    section=4;
end

if OverArclength || strcmp(Quantity, 'Cfint')
    % arclength
    x= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(wake.s)) + wake.s];
    xstr='s';
else
    x= [prf.nodes.X'; wake.x];
    xstr='x';
end

% chose section
if section==1
    ind=(prf.M:-1:1);              % suction side node indizes 
elseif section==2
    ind=(prf.M:prf.N);             % pressure side node indizes 
elseif section==3
    ind=(prf.N + 1:prf.N + wake.N);    % wake node indizes
elseif section==4
    ind1=(prf.M:-1:1); 
    ind2=(prf.M:prf.N);  
else
    ind1=(prf.M:-1:1); 
    ind2=(prf.M:prf.N + wake.N);  
end

%legend for sections
if section==4
    legendstr={'suction side','pressure side'}; 
elseif section==5
    legendstr={'suction side','pressure side / wake'}; 
end

% plots with "special" threatment: Cfint, delta
plotted=false;
if strcmp(Quantity ,'Cfint') || strcmp(Quantity ,'cfint')
    %------------------------- Cfint -----------------------------------------------
    figure
    hold on
    title('integral wall friction coefficient')
    xlabel(xstr)
    ylabel( ' [ C_f] ' )
    
    sU=[0,prf.sU(end:-1:1)];
    sL=[0,prf.sL];
    if section==1
        plot( sU, sol.CI_U );
    elseif section==2
        plot( sL, sol.CI_L );
    else
        plot( sU, sol.CI_U );
        plot( sL, sol.CI_L );
        legend('suction side','pressure side','location','best'); 
    end
    plotted=true;    
elseif strcmp(Quantity ,'delta') 
    %------------------------- delta  -----------------------------------------------
    figure
    hold on
    title('Displacement thickness and momentum thickness')
    xlabel(xstr)
    ylabel( ' \delta ' )
    
    if section<4
        plot( x(ind), sol.D(ind) );
        plot( x(ind), sol.T(ind) );
        legend('\delta_1','\delta_2','location','best'); 
    else
        plot( x(ind1), sol.D(ind1) ,'b');
        plot( x(ind2), sol.D(ind2) ,'r');
        plot( x(ind1), sol.T(ind1) ,'b --');
        plot( x(ind2), sol.T(ind2) ,'r --');
        legend('$\delta_1$ suction side','$\delta_1$ pressure side','$\delta_2$ suction side','$\delta_2$ pressure side','location','best'); 
    end  
    plotted=true;
end


% determine title, y-axis label and y quantity
if strcmp(Quantity ,'tau') || strcmp(Quantity ,'Tau')
    titlestr='Wall shear stress';
    ystr=' \tau_w / \rho U_\infty^2 ';
    y_arr=sol.tau;
elseif strcmp(Quantity ,'Cf') || strcmp(Quantity ,'cf')
    titlestr='wall friction coefficient normalized with U';
    ystr=' C_f=2 \tau_w / \rho U^2 ';
    y_arr=sol.cf;
elseif strcmp(Quantity ,'Cp') || strcmp(Quantity ,'cp')
    titlestr='pressure coefficient';
    ystr=' C_p=1 - U^2 /U_\infty^2 ';
    y_arr=sol.Cp;    
elseif strcmp(Quantity ,'U') || strcmp(Quantity ,'u')    
    titlestr='boundary layer edge velocity';
    ystr='U';
    y_arr=sol.U;
elseif strcmp(Quantity ,'CD') || strcmp(Quantity ,'cD')   
    titlestr='Dissipation coefficient';
    ystr= 'C_D';
    y_arr=  sol.CD;  
elseif strcmp(Quantity ,'D') 
    titlestr='Dissipation integral';
    ystr=' D ';
    y_arr=  sol.CD.*sol.U.^3 ;    
elseif strcmp(Quantity ,'H12')     
    titlestr='Shapeparameter H12';
    ystr=' H_1_2 ' ;
    y_arr=  sol.HK; 
elseif strcmp(Quantity ,'H32')     
    titlestr='Shapeparameter H32';
    ystr=' H_3_2 ' ;
    y_arr=  sol.HS;     
elseif strcmp(Quantity ,'Ret') || strcmp(Quantity ,'ReT')
    titlestr='Shapeparameter H32';
    ystr=' Re_T ' ;
    y_arr= sol.Ret;     
elseif strcmp(Quantity ,'q') || strcmp(Quantity ,'Q')    
    titlestr='source contribution';
    ystr=' q ' ;
    q=[FiniteDifferences(sol.m(prf.Nle-1:-1:1),prf.sU(end:-1:1)) ;...
       FiniteDifferences(sol.m(prf.Nle:prf.N),prf.sL)];
    y_arr= q; 
elseif strcmp(Quantity ,'m') || strcmp(Quantity ,'M')    
    titlestr='mass defect';
    ystr=' m ' ;
    y_arr= sol.m; 
elseif strcmp(Quantity ,'n') ||strcmp(Quantity ,'N')
    titlestr='amplification exponent';
    ystr=' n ' ;
    numUP = prf.M-sol.iTran(1)+1;
    numLO = sol.iTran(2)-prf.M;
    if section==1
        ind=ind(1:numUP);
    elseif section==2
        ind=ind(2:numLO);
    elseif section==3
        disp('no amplification in wake')
        return
    elseif section==4 || section==5
        ind1=ind1(1:numUP);
        ind2=ind2(2:numLO);
    end
    n=sol.c;
    n(sol.iTran(1))=sol.tran.n2(1);
    n(sol.iTran(2))=sol.tran.n2(2);
    n(1:sol.iTran(1)-1)=0;
    n(sol.iTran(2):end)=0;
    y_arr= n; 
    
    nkrit= (sol.tran.n2(1)-n(sol.iTran(1)+1))*sol.tran.Llam(1)...
        /(sol.tran.Llam(1) + sol.tran.Lturb(1))+n(sol.iTran(1)+1);  
elseif strcmp(Quantity ,'c') ||strcmp(Quantity ,'C')||strcmp(Quantity ,'Ctau')||strcmp(Quantity ,'ctau')
    titlestr='shear stress coefficient';
    ystr='C_\tau' ;
    c=sol.c;
    c(sol.iTran(1)+1:sol.iTran(2)-1)=0;
    y_arr= c; 
end
    

% Do the plots
    
if ~plotted
   figure 
   hold on
   title(titlestr) 
   xlabel(xstr)
   ylabel(ystr)
    if section<4
        plot( x(ind), y_arr(ind) );
    else
        plot( x(ind1), y_arr(ind1) );
        plot( x(ind2), y_arr(ind2) );
    end
   if section>3     
        legend(legendstr,'location','best');  
   end
end


% set limits
if ~(strcmp(Quantity ,'Cfint') || strcmp(Quantity ,'cfint'))
    if section<4
        xmin=min(x(ind));
        xmax=max(x(ind));
        
    else
        xmin=min(x(ind1));
        xmax=max(x(ind2));
    end
    xlim([xmin,xmax])
end


if strcmp(Quantity ,'n') ||strcmp(Quantity ,'N')
   line([xmin  xmax]  , [nkrit  nkrit],'color','black');%,'LineStyle','--'); 
   text(0.5*(xmin+xmax), nkrit*0.96, 'n_k_r_i_t');
end

%-------------------------------------------------------------------------------------------------------------------------

% % plot "thicked" shape
% %---------------------------------------------------
% XT=prf.nodes.X - prf.nodes.n(1,:).*transpose(sol.D(1:prf.N));
% YT=prf.nodes.Y - prf.nodes.n(2,:).*transpose(sol.D(1:prf.N));
% XTw=wake.x' + wake.nn(1,:).*transpose(sol.D(prf.N+1:end));
% 
% %weightening for Wake 
% wUpper= (0.5* prf.gap + sol.D(1))/sol.D(prf.N+1);
% wLower= (0.5* prf.gap + sol.D(prf.N))/sol.D(prf.N+1);
% fU=1+wUpper -wLower;
% fL=1-wUpper +wLower;
% 
% % makes continous transition onto wake
% nwU= [-prf.nodes.n(:,1), wake.nn(:,2:end)];
% nwL= [prf.nodes.n(:,end), wake.nn(:,2:end)];
% XUw=[XT(1),XTw(2:end)];
% XLw=[XT(end),XTw(2:end)];
% 
% YTw =wake.y' + fU* nwU(2,:).*transpose(sol.D(prf.N+1:end))/2;
% YTw2=wake.y' - fL* nwL(2,:).*transpose(sol.D(prf.N+1:end))/2;
% 
% xtt= [XUw(end:-1:1),XT,XLw];
% ytt= [YTw(end:-1:1),YT,YTw2];
% 
% tmp=round(prf.alfa*180/pi);
% str={['\alpha=',num2str(tmp),' degree'],['Re=',num2str(Re)],['C_L=',num2str(sol.CL)],['C_d=',num2str(sol.Cdrag)],['C_\nu=',num2str(sol.Cnu)],['C_L/C_d=',num2str(sol.CL/sol.Cdrag)] };
% 
% figure; 
% hold on; box on;
% plot([prf.panels.X]',[prf.panels.Y]','k','Linewidth',2);
% plot(wake.x,wake.y,'k','Linewidth',1.2);
% plot(xtt,ytt,'b');
% % plot(XUw,YTw,'b');
% % plot(XLw,YTw2,'b');
% text(wake.x(end)-0.7*(wake.x(end)-prf.nodes.X(1)) ,0.48, str);
% title('Profile with displacement thickness')
% axis equal; xlabel('x'); ylabel('y') 
% clear XT YT XTw wUpper wLower fU fL nwU nwL XUw XLw YTw YTw2
% %---------------------------------------------------



end



% dat=[xtt;ytt]';
% 
% exp=dat(1:3:end,:);
% dlmwrite('thicked',exp);






