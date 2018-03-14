function [wake] = GetWakeStreamline( fld,prf,NW )
%GETWAKESTREAMLINE  calculates the nodes on the wake streamline. 
%                   uses Newton method to find points with psi=psi0
%                   NW: number of wake nodes



%        calculate first wake tangent vector
%           -> mean of tangent vectors of first and N-th panel 
e1=prf.panels.e(:,1);
eN=prf.panels.e(:,end-1);
s=0.5*transpose(eN-e1);
%if s(1)<0; s=-s; end;



if prf.IsSharp
    xTE=[prf.panels.X(2,end) prf.panels.Y(2,end)]+0.00001*s;
    
    psi0=evaluateInviscFieldSol(xTE,fld,prf);
else
    %    calculate TE midpoint as first wake node
    x1=[prf.panels.X(1,1) prf.panels.Y(1,1)]; 
    xN=[prf.panels.X(1,end) prf.panels.Y(1,end)];

    xTE= (x1+xN)/2 + 0.0001*s; % Midpoint of TE
    psi0=fld.psi0;
    %psi0=evaluateInviscFieldSol(xTE,fld,prf);
end

% starting step size of wake equals the Panel length of first panel
h=0.5*( prf.panels.L(1) + prf.panels.L(prf.N-1) );

grading=Grading(h,prf.c,NW);


% guess second wake point 
guess= xTE + h*s;

xw=zeros(NW,2); xw(1,:)=xTE;
sw=zeros(1,NW);
lw=zeros(1,NW-1);


for i=1:NW-1 
    xw(i+1,:)=guess;
  
    %------------- iterate to find Loadingpoint with exact psi=psi0  -------
    res=1;k=0; dpmin=2;noConv=false;
    while res>1e-6 && k<40 % prevent endless loop
        psi= evaluateInviscFieldSol(xw(i+1,:),fld,prf);
        
        
        [dg_dx,dg_dy]=GradPsi(xw(i+1,1), xw(i+1,2),prf);
        dpsi_dx=-dg_dx*fld.gamma;
        dpsi_dy=-dg_dy*fld.gamma;

        ht=norm(xw(i+1,:)-xw(i,:));
        % Newton-Method

        J=[dpsi_dx, dpsi_dy; 2*(xw(i+1,1)-xw(i,1)), 2*(xw(i+1,2)-xw(i,2))];
        f=[psi-psi0, ht^2-h^2];

        dx=  (f(2)*J(1,2) - f(1)*J(2,2))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );
        dy=  (f(1)*J(2,1) - f(2)*J(1,1))/( J(1,1)*J(2,2)- J(1,2)*J(2,1) );
        
        
        
        % relativ difference
        dp= abs((psi-psi0)/psi0);
        if dp<dpmin
           dpmin=dp; % minimum deviation of psi
           xmin=xw(i+1,:);
        end
        
        res=max(abs(f./[psi0,h]));
        % refresh for new iteration
        xw(i+1,:)= xw(i+1,:) + [dx, dy];
        
        % if Newtonmethod does not converge or converges to wrong solution
        if xw(i+1,1)<xw(i,1) || abs(xw(i+1,2)-xw(i,2)) > 0.5*h
           noConv=true; break 
        end
        
        k=k+1;
    end

    if noConv || dpmin>0.008 
        % in case of divergence or to big deviation of psi0 -> try intersection method
        [xI, psiI]=Intersection(fld,prf,guess,h,psi0);
        if abs((psiI-psi0)/psi0)<0.008 ;
           xw(i+1,:)= xI; 
        else
           xw(i+1,:)=guess; 
        end
    else
        xw(i+1,:)= xmin;
    end
    
    lw(i)=norm(xw(i+1,:)-xw(i,:));
    delta=xw(i+1,:)-xw(i,:);
    guess= xw(i+1,:) + grading*delta;
    sw(i+1)=sw(i) + lw(i);
    h=lw(i)*grading;
    
    clear nearest
end


%tangential vectors

ew=transpose(xw(2:end,:)-xw(1:end-1,:)); ew=[ew(1,:)./lw; ew(2,:)./lw];
nw=[-ew(2,:);ew(1,:)];

% normal vector of node -> mean value of panel normal vectors
nwn=(nw(:,1:end-1)+nw(:,2:end))/2;
nwn=[nw(:,1),nwn, nw(:,end)];

wake.theta=atan2(ew(1,:),-ew(2,:));
wake.x=xw(:,1);
wake.y=xw(:,2);
wake.e=ew;
wake.n=nw;
wake.L=lw;
wake.s=sw;
wake.nn=nwn;
wake.N=NW;


% Trailing edge Gap

AN=prf.ScrossT.*prf.panels.L(end);
ZN=1-wake.s/(2.5*AN);
xp1=prf.nodes.X(2)-prf.nodes.X(1);
xpN=prf.nodes.X(end)-prf.nodes.X(end-1);
yp1=prf.nodes.Y(2)-prf.nodes.Y(1);
ypN=prf.nodes.Y(end)-prf.nodes.Y(end-1);
Cross= (xp1*ypN-yp1*xpN)/sqrt( (xp1^2+yp1^2)*(xpN^2+ypN^2) );
D=Cross/sqrt(1-Cross^2);
D=max(D,-3/2.5);
D=min(D, 3/2.5);
A= 3+2.5*D;
B=-2-2.5*D;
wg=zeros(size(wake.s));
wg(ZN>0)= AN* (A+B*ZN(ZN>0)).*ZN(ZN>0).^2;
wake.gap=wg';


% % starting y-step size
% dy1=0.2*( prf.nodes.Y(2)-prf.nodes.Y(end-1) );
% for i=1:NW-1 
%     psi1= evaluateInviscFieldSol(guess,fld,prf);
%     
%     
%     psi2= evaluateInviscFieldSol(guess + [0, dy1],fld,prf);
%     
%     dy= (psi0-psi1)/(psi2-psi1)*dy1;
%     xw(i+1,:)=guess;
%     
%     res=1;k=0;
%     while res>1e-9 && k<20 % prevent endless loop
%         psi2= evaluateInviscFieldSol(xw(i+1,:) + [0, dy],fld,prf);
%         dy= (psi0-psi1)/(psi2-psi1)*dy;
%         
%         % refresh for new iteration
%         xw(i+1,2)= xw(i+1,2) + dy;
%         psi1=psi2;
%         k=k+1;
%     end
%     
%     
%     %grading=1+(grading-1)/(1+0.0002*(i-1)^2); % make grading go to 1 at end of wake
%     guess=xw(i+1,:) + grading*(xw(i+1,:)-xw(i,:));
%     lw(i)=norm(xw(i+1,:)-xw(i,:));
%     sw(i+1)=sw(i) + lw(i);
% end
end

