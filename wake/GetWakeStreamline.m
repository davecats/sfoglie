function [wake] = GetWakeStreamline( fld,prf,NW ,grading)
%GETWAKESTREAMLINE  calculates the noedes on the wake streamline. 
%                   Finds points where psi-psi0 goes to zero -> Sekantenmethode
%                   NW: number of wake nodes
%                   grading: default set to 1.1


if nargin==3
grading=1.1;
end


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

    xTE= (x1+xN)/2; % Midpoint of TE
    psi0=fld.psi0;
    %psi0=evaluateInviscFieldSol(xTE+0.0001*s,fld,prf);
end




% starting step size of wake equals the Panel length of first panel
h=prf.panels.L(1);

% guess second wake point 
guess= xTE+1.01*h*s;

xw=zeros(NW,2); xw(1,:)=xTE;
sw=zeros(1,NW);
lw=zeros(1,NW-1);
%value of streamfunktion for guessing point

dy=30*abs(h*s(2));

for i=1:NW-1 % each node of wake
    psig=evaluateInviscFieldSol(guess,fld,prf);
    found=false;k=1;
    
    % find point on wake streamline -> psi=psi0
    
    % Einschlussintervall finden
    %------------------------------------------------
    while found==false && k<20 % prevent endless loop
        psip = evaluateInviscFieldSol(guess+[0, k*dy],fld,prf);
        d1=psig-psi0;
        d2=psip-psi0;
        if sign(d1)*d2 < 0 % Einschlussintervall gefunden
          psi1=psig;
          psi2=psip;
          found=true;
          y2 = guess(2) + k*dy;
        else
          psim=evaluateInviscFieldSol(guess-[0, k*dy],fld,prf);  
          d2=psim-psi0;  
          if sign(d1)*d2 < 0 % Einschlussintervall gefunden
            psi1=psig;
            psi2=psim;  
            found=true;
            y2 = guess(2) - k*dy;
          end
        end
        k=k+1;   
    end
    %------------------------------------------------
    y1=guess(2);
    dn=1; l=0;
    % Sekanten Methode
    %------------------------------------------------
    res=1e-12; % residuum
    while  abs(dn)>res && l<40 % prevent endless loop
        d1=psi1-psi0;
        d2=psi2-psi0;
        yn=(abs(d1)*y2 + abs(d2)*y1) /(abs(d1)+abs(d2));
        new=[guess(1), yn];
        psin=evaluateInviscFieldSol(new,fld,prf);
        dn=psin-psi0;
        if sign(dn)*d1 <0
           psi2 = psin;
           y2=yn;
        else
           psi1 = psin;
           y1=yn; 
        end
        l=l+1;
    end
    %------------------------------------------------
    xw(i+1,:)=[guess(1), yn];
    grad=1+(grading-1)/(1+0.0002*(i-1)^2); % make grading go to 1 at end of wake
    guess=xw(i+1,:) + grad*(xw(i+1,:)-xw(i,:));
    lw(i)=norm(xw(i+1,:)-xw(i,:));
    sw(i+1)=sw(i) + lw(i);
    
%     if lw(i)>20*lw(1)
%        grading=1; 
%     end
    
end

%tangential vektors

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
wake.GAP=wg;


end

