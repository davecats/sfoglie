function [fig] = PlotProfile( prf,wake,sol, mode) 
%PLOTPROFIL  mode=1: Plots the profile and wake with thicked shape due to displacement
%            mode=2: Plots the profile and wake with node positions
%            mode=3: Plots Profile with Transitionpoints and blowing distribution   

if nargin==3
    mode=1;
end


% plot "thicked" shape
%---------------------------------------------------
if mode==1
    Re=evalin('base','Re');
    
    
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

    fig=figure; 
    hold on; box on;
    plot([prf.panels.X],[prf.panels.Y],'k','Linewidth',1.5);
    plot(wake.x,wake.y,'k','Linewidth',1.2);
    plot(xtt,ytt,'b');
    % plot(XUw,YTw,'b');
    % plot(XLw,YTw2,'b');
    text(wake.x(end)-0.7*(wake.x(end)-prf.nodes.X(1)) ,0.48*prf.c, str);
    title('Profile with displacement thickness')
    axis equal; xlabel('x'); ylabel('y') 
    clear XT YT XTw wUpper wLower fU fL nwU nwL XUw XLw YTw YTw2
    %---------------------------------------------------
    % dat=[xtt;ytt]';
    % 
    % exp=dat(1:3:end,:);
    % dlmwrite('thicked',exp);
elseif  mode ==2
    % plot nodes of airfoil geometrie
    
    x= [wake.x(1); prf.nodes.X';wake.x];
    y= [wake.y(1); prf.nodes.Y';wake.y];
    
    nn=[prf.nodes.n,wake.nn(:,2:end)]';
    lls=0.01*prf.c;
    
    PP=[[prf.nodes.X';wake.x(2:end)] , [prf.nodes.Y';wake.y(2:end) ]];
    startP= PP + lls*nn;
    endP  = PP - lls*nn;
    
    
    str={['N=',num2str(prf.N)],['N_w=',num2str(wake.N)]};
    
    fig=figure; 
    hold on; box on;
    plot(x,y,'k','Linewidth',1);
    for i=1:length(PP)
        line([startP(i,1) endP(i,1)]  , [startP(i,2) endP(i,2)],'color','black','Linewidth',0.7);
    end
    text(wake.x(end)-0.7*(wake.x(end)-prf.nodes.X(1)) ,0.6*prf.c, str);
    title('Profile with nodes')
    axis equal; xlabel('x'); ylabel('y') 

elseif mode==3 || mode==4
    
    xUtr= [prf.nodes.X(sol.iTran(1));prf.nodes.Y(sol.iTran(1))] + sol.tran.Lturb(1)*prf.panels.e(:,sol.iTran(1));
    xLtr= [prf.nodes.X(sol.iTran(2));prf.nodes.Y(sol.iTran(2))] - sol.tran.Lturb(2)*prf.panels.e(:,sol.iTran(2)-1);
    lls=0.03*prf.c;
    endU= xUtr - lls*prf.panels.n(:,sol.iTran(1));
    endL= xLtr - lls*prf.panels.n(:,sol.iTran(2)-1);
    
    str={'Transition'...
         ['suction side  : x/c=',num2str(round(xUtr(1),3)), '  y/c=',num2str(round(xUtr(2),3))]...
         ['pressure side: x/c=',num2str(round(xLtr(1),3)), '  y/c=',num2str(round(xLtr(2),3))]};
       
    fig=figure; 
    hold on; box on;
    l1=line([xUtr(1) endU(1)]  , [xUtr(2) endU(2)],'color','r','Linewidth',0.7);
    line([xLtr(1) endL(1)]  , [xLtr(2) endL(2)],'color','r','Linewidth',0.7);

    ind=find(abs(sol.Vb(1:prf.N) )>1e-7); 
    if ~isempty(ind)
        add=' and blowing profiles';
        Faktor= 0.02*prf.c/max(sol.Vb(1:prf.N));
        vb= [prf.nodes.X;prf.nodes.Y] - Faktor*([1;1]*sol.Vb(1:prf.N)' ).*prf.nodes.n;
        l2=plot(vb(1,:),vb(2,:),'b');
        legendE1='TP';
        if mode==4; legendE1=[legendE1,' blowing case']; end
        legend([l1 l2],{legendE1,'blowing velocity'},'location','northwest')
        vTT=vb(:,ind);
        for i=1:length(vTT(1,:))
            line([prf.nodes.X(ind(i)) vTT(1,i)]  , [prf.nodes.Y(ind(i)) vTT(2,i)],'color','b');
        end
    else
        legend(l1,{'transition point'},'location','northwest')
        add='';
    end  
    plot([prf.panels.X],[prf.panels.Y],'k','Linewidth',1)

    if mode==3; text(0.4*prf.c, 0.3*prf.c, str); end
    title(['Profile with Transition',add])
    axis equal; xlabel('x'); ylabel('y')  
     
end






end




