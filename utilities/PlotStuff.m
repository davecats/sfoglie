function [  ] = PlotStuff( prf,wake,sol ,Val ,ind )
%PLOTSTUFF Function for Plotting 


if nargin== 4
   ind=1:prf.N+wake.N; 
end

tr=false;
Re=evalin('base','Re');
nkrit=evalin('base','nkrit');

% if Re== 400000
%     str1= 'Re4e5_';
% elseif abs(Re-65360)<1
%     str1= 'Re65360_';
% elseif  abs(Re-100000)<1
%     str1= 'Re1e5_';
% end
% tmp=round(prf.alfa*180/pi);
% if abs(tmp-2)<2e-3
%     str2='alfa2_';
% elseif abs(tmp-5)<2e-3
%     str2='alfa5_';tr=true;
% %elseif abs(tmp)<2e-3
% else
%     str2='alfa0_';   
% end
% 
% 
% if abs(nkrit-0.105)<2e-3
%     str3='N0105';
% elseif abs(nkrit-0.155)<2e-3
%     str3='N0155';
% elseif abs(nkrit-0.15)<2e-3
%     str3='N015';   
% %elseif abs(nkrit-9)<2e-3
% else
%     str3='N9'; 
% end
% if tr; str3='Trip01'; end
% str=[str1,str2,str3];

% % s, x, y, U, D, T, Cf, H
% data =importdata(['./XFoilWerte/',str]);
% data=data.data;
% 
% sX=data(:,1);
% XX=data(:,2);
% 
% UX=abs(data(:,4));
% DX=data(:,5);
% TX=data(:,6);
% tauX=data(:,7)*0.5;
% HX=DX./TX;
% ReX=UX.*TX/nu;

% 
% NX=length(XX<=1);
% NleX= find(UX<min(UX)+1e-4);
% 
% CfX=2*tauX./UX.^2;

%sges= [prf.nodes.X'; wake.x];
%prf.nodes.X

% data=load('./XFoilWerte/MRDU.txt');
% data=[data(1:2:end-1,:), data(2:2:end,:)];
% data(1:prf.Nle-1,:)=data(prf.Nle-1:-1:1,:);
% UX=abs(data(:,4));
% DX=data(:,3);
% TX=data(:,2);
% tauX=data(:,4);

sges= [prf.s, (prf.s(end)+prf.panels.L(end)/2)*ones(size(wake.s)) + wake.s]; 

if strcmp(Val,'delta')
    figure
    hold on
    plot( sges(ind), sol.D(ind) );
    plot( sges(ind), sol.T(ind) );
    %plot( sX(ind), DX(ind) );
    %plot( sX(ind), TX(ind) );
    %legend('$\delta_1$','$\delta_2$','$\delta_1$ XFoil','$\delta_2$ XFoil','location','best'); 
    title('Displacement thickness and momentum thickness')
    xlabel(' s')
    ylabel( '\delta_i' )
    legend('$\delta_1$','$\delta_2$','location','best'); 
elseif strcmp(Val,'U')
    figure
    hold on
    plot( sges(ind), sol.U(ind) );
%     plot( sges(ind), UX(ind) );
    title('Boundary layer edge velocity')
    xlabel(' s')
    ylabel( ' U=u|_\delta ' )
    legend('Lsg','XFoil','location','best'); 
elseif strcmp(Val,'cf') || strcmp(Val,'Cf')
    figure
    hold on
    plot( prf.nodes.X(1:prf.Nle), sol.Cf(1:prf.Nle));
    plot( prf.nodes.X(prf.Nle:prf.N), sol.Cf(prf.Nle:prf.N));
%     plot( XX(1:NleX-1), CfX(1:NleX-1));
%     plot( XX(NleX-1:NX), CfX(NleX-1:NX));
%     legend('suction side','pressure side','XFoil: suction','XFoil: pressure','location','best');
    title('skin friction coefficient')
    xlabel(' x ')
    ylabel( ' C_f ' )
    legend('suction side','pressure side','location','best');
elseif strcmp(Val,'tau') || strcmp(Val,'Tau')
    figure
    hold on
    plot( prf.nodes.X(1:prf.Nle), sol.tau(1:prf.Nle));
    plot( prf.nodes.X(prf.Nle:prf.N), sol.tau(prf.Nle:prf.N));
%     plot( XX(1:NleX-1), tauX(1:NleX-1));
%     plot( XX(NleX:NX), tauX(NleX:NX));
%     legend('suction side','pressure side','XFoil: suction','XFoil: pressure','location','best'); 
    title('Wall shear stress')
    xlabel(' x ')
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
    plot( prf.nodes.X(1:prf.Nle), sol.Cp(1:prf.Nle));
    plot( prf.nodes.X(prf.Nle:prf.N), sol.Cp(prf.Nle:prf.N));
    legend('suction side','pressure side','location','best')
    title('Pressure coefficient')
    xlabel('x')
    ylabel( ' C_p ' )    
end    
    
%-------------------------------------------------------------------------------------------------------------------------

% plot "thicked" shape
%---------------------------------------------------
XT=prf.nodes.X - prf.nodes.n(1,:).*transpose(sol.D(1:prf.N));
YT=prf.nodes.Y - prf.nodes.n(2,:).*transpose(sol.D(1:prf.N));
XTw=wake.x' + wake.nn(1,:).*transpose(sol.D(prf.N+1:end));

%weightening for Wake 
wUpper= (0.5* prf.gap + sol.D(1))/sol.D(prf.N+1);
wLower= (0.5* prf.gap + sol.D(prf.N))/sol.D(prf.N+1);
fU=1+wUpper -wLower;
fL=1-wUpper +wLower;

% makes continous transition onto wake
nwU= [-prf.nodes.n(:,1), wake.nn(:,2:end)];
nwL= [prf.nodes.n(:,end), wake.nn(:,2:end)];
XUw=[XT(1),XTw(2:end)];
XLw=[XT(end),XTw(2:end)];

YTw =wake.y' + fU* nwU(2,:).*transpose(sol.D(prf.N+1:end))/2;
YTw2=wake.y' - fL* nwL(2,:).*transpose(sol.D(prf.N+1:end))/2;

xtt= [XUw(end:-1:1),XT,XLw];
ytt= [YTw(end:-1:1),YT,YTw2];

tmp=round(prf.alfa*180/pi);
str={['\alpha=',num2str(tmp),' degree'],['Re=',num2str(Re)],['C_L=',num2str(sol.CL)],['C_d=',num2str(sol.Cdrag)],['C_\nu=',num2str(sol.Cnu)],['C_L/C_d=',num2str(sol.CL/sol.Cdrag)] };

figure; 
hold on; box on;
plot([prf.panels.X]',[prf.panels.Y]','k','Linewidth',2);
plot(wake.x,wake.y,'k','Linewidth',1.2);
plot(xtt,ytt,'b');
% plot(XUw,YTw,'b');
% plot(XLw,YTw2,'b');
text(wake.x(end)-0.7*(wake.x(end)-prf.nodes.X(1)) ,0.48, str);
title('Profile with displacement thickness')
axis equal; xlabel('x'); ylabel('y') 
clear XT YT XTw wUpper wLower fU fL nwU nwL XUw XLw YTw YTw2
%---------------------------------------------------



end



% dat=[xtt;ytt]';
% 
% exp=dat(1:3:end,:);
% dlmwrite('thicked',exp);






