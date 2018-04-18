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

if strcmp(Quantity ,'tau') || strcmp(Quantity ,'Tau')
    %-------------------------  TAU -----------------------------------------------
    figure
    hold on
    title('Wall shear stress')
    xlabel(xstr)
    ylabel( ' \tau_w / \rho U_\infty^2 ' )
    
    if section<4
        plot( x(ind), sol.tau(ind) );
    else
        plot( x(ind1), sol.tau(ind1) );
        plot( x(ind2), sol.tau(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end
    
    
elseif strcmp(Quantity ,'Cfint') || strcmp(Quantity ,'cfint')
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

    
elseif strcmp(Quantity ,'Cf') || strcmp(Quantity ,'cf')
    %------------------------- Cf   -----------------------------------------------
    figure
    hold on
    title('wall friction coefficient normalized with U')
    xlabel(xstr)
    ylabel( ' C_f=2 \tau_w / \rho U^2 ' )
    
    if section<4
        plot( x(ind), sol.Cf(ind) );
    else
        plot( x(ind1), sol.Cf(ind1) );
        plot( x(ind2), sol.Cf(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end
    
elseif strcmp(Quantity ,'Cp') || strcmp(Quantity ,'cp')
    %------------------------- Cp   -----------------------------------------------
    figure
    hold on
    title('pressure coefficient')
    xlabel(xstr)
    ylabel( ' C_p=1 - U^2 /U_\infty^2 ' )
    
    if section<4
        plot( x(ind), sol.Cp(ind) );
    else
        plot( x(ind1), sol.Cp(ind1) );
        plot( x(ind2), sol.Cp(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end
    
    
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
    
    
elseif strcmp(Quantity ,'U') || strcmp(Quantity ,'u')
    %------------------------- U  -----------------------------------------------
    figure
    hold on
    title('boundary layer edge velocity')
    xlabel(xstr)
    ylabel( ' U ' )
    
    if section<4
        plot( x(ind), sol.U(ind) );
    else
        plot( x(ind1), sol.U(ind1) );
        plot( x(ind2), sol.U(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end  

    
elseif strcmp(Quantity ,'CD') || strcmp(Quantity ,'cD')
    %------------------------- CD  -----------------------------------------------
    figure
    hold on
    title('Dissipation coefficient')
    xlabel(xstr)
    ylabel( ' C_D ' )
    
    if section<4
        plot( x(ind), sol.CD(ind) );
    else
        plot( x(ind1), sol.CD(ind1) );
        plot( x(ind2), sol.CD(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end     
    
elseif strcmp(Quantity ,'D') 
    %------------------------- D  -----------------------------------------------
    figure
    hold on
    title('Dissipation integral')
    xlabel(xstr)
    ylabel( ' D ' )
    
    if section<4
        plot( x(ind), sol.CD(ind).*sol.U(ind).^3 );
    else
        plot( x(ind1), sol.CD(ind1).*sol.U(ind1).^3 );
        plot( x(ind2), sol.CD(ind2).*sol.U(ind2).^3 );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end     
    
elseif strcmp(Quantity ,'H12') 
    %------------------------- H12  -----------------------------------------------
    figure
    hold on
    title('Shapeparameter H12')
    xlabel(xstr)
    ylabel( ' H_1_2 ' )
    
    if section<4
        plot( x(ind), sol.HK(ind) );
    else
        plot( x(ind1), sol.HK(ind1) );
        plot( x(ind2), sol.HK(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end      
    
    
elseif strcmp(Quantity ,'H32') 
    %------------------------- H32  -----------------------------------------------
    figure
    hold on
    title('Shapeparameter H32')
    xlabel(xstr)
    ylabel( ' H_3_2 ' )
    
    if section<4
        plot( x(ind), sol.HS(ind) );
    else
        plot( x(ind1), sol.HS(ind1) );
        plot( x(ind2), sol.HS(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end        
    
 elseif strcmp(Quantity ,'Ret') || strcmp(Quantity ,'ReT')
    %------------------------- Ret  -----------------------------------------------
    figure
    hold on
    title('momentumthickness based Reynolds-number')
    xlabel(xstr)
    ylabel( ' Re_T ' )
    
    if section<4
        plot( x(ind), sol.Ret(ind) );
    else
        plot( x(ind1), sol.Ret(ind1) );
        plot( x(ind2), sol.Ret(ind2) );
        if section==4
            legend('suction side','pressure side','location','best'); 
        else
             legend('suction side','pressure side / wake','location','best'); 
        end
    end     
    
end

% limits
if section<4
    xlim([min(x(ind)),max(x(ind))])
else
    xlim([min(x(ind2)),max(x(ind2))])
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






