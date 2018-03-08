function [  ] = PlotStuff( prf,wake,sol ,Val ,ind )
%PLOTSTUFF Function for Plotting 


if nargin== 4
   ind=1:prf.N+wake.N; 
end

tr=false;
Re=evalin('base','Re');
nkrit=evalin('base','nkrit');


sges= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(wake.s)) + wake.s]; 



if strcmp(Val,'delta')
    figure
    hold on
    plot( sges(ind), sol.D(ind) );
    plot( sges(ind), sol.T(ind) );
    title('Displacement thickness and momentum thickness')
    xlabel(' s')
    ylabel( '\delta_i' )
    legend('$\delta_1$','$\delta_2$','location','best'); 
elseif strcmp(Val,'U')
    figure
    hold on
    plot( sges(ind), sol.U(ind) );
    title('Boundary layer edge velocity')
    xlabel(' s')
    ylabel( ' U=u|_\delta ' )
    legend('Lsg','location','best'); 
elseif strcmp(Val,'cf') || strcmp(Val,'Cf')
    figure
    hold on
    plot( prf.nodes.X(1:prf.Nle-1), sol.Cf(1:prf.Nle-1));
    plot( prf.nodes.X(prf.Nle:prf.N), sol.Cf(prf.Nle:prf.N));
    title('skin friction coefficient')
    xlabel(' s')
    ylabel( ' C_f ' )
    legend('suction side','pressure side','location','best'); 
elseif strcmp(Val,'tau') || strcmp(Val,'Tau')
    figure
    hold on
    plot( prf.nodes.X(1:prf.Nle-1), sol.tau(1:prf.Nle-1));
    plot( prf.nodes.X(prf.Nle:prf.N), sol.tau(prf.Nle:prf.N));
    title('Wall shear stress')
    xlabel(' s')
    ylabel( ' \tau_w ' )
    legend('suction side','pressure side','location','best'); 
elseif strcmp(Val,'CD') || strcmp(Val,'cD')
    figure
    hold on
    plot( sges(ind), sol.CD(ind) );
    title('Dissipation coefficient')
    xlabel(' s')
    ylabel( ' C_D ' )
 elseif strcmp(Val,'D') 
    figure
    hold on
    plot( sges(ind), sol.CD(ind).*sol.U.^2 );
    title('Dissipation integral')
    xlabel(' s')
    ylabel( ' D ' )   
  elseif strcmp(Val,'cp') || strcmp(Val,'Cp')
    figure
    hold on
    plot( prf.nodes.X(1:prf.Nle-1), sol.Cp(1:prf.Nle-1));
    plot( prf.nodes.X(prf.Nle:prf.N), sol.Cp(prf.Nle:prf.N));
    legend('suction side','pressure side','location','best')
    title('Pressure coefficient')
    xlabel('x')
    ylabel( ' C_p ' )    
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
% tmp=round(prf.alfa*180/pi);
% str={['\alpha=',num2str(tmp),' degree'],['Re=',num2str(Re)],['C_L=',num2str(sol.CL)] };
% 
% figure; 
% hold on; box on;
% plot([prf.panels.X]',[prf.panels.Y]','k','Linewidth',2);
% plot(wake.x,wake.y,'k','Linewidth',1.2);
% plot(XT,YT,'b');
% plot(XUw,YTw,'b');
% plot(XLw,YTw2,'b');
% text(wake.x(end)-0.5*(wake.x(end)-prf.nodes.X(1)) ,0.48, str);
% title('Profile with displacement thickness')
% axis equal; xlabel('x'); ylabel('y') 
% clear XT YT XTw wUpper wLower fU fL nwU nwL XUw XLw YTw YTw2
% %---------------------------------------------------



end

