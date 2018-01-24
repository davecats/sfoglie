function LEstr = LEstreamline( fld,prf,N ,grading)
%LESTREAMLINE Summary of this function goes here
%   Detailed explanation goes here



if nargin==3
grading=1.15;
end


s=-transpose( prf.panels.n(:,prf.Nle-1));
xLE=[prf.panels.X(1,prf.Nle-1), prf.panels.Y(1,prf.Nle-1)] + prf.LE1*transpose(prf.panels.e(:,prf.Nle-1)); 

psi0=evaluateInviscFieldSol(xLE,fld,prf);


% starting step size of wake equals the Panel length of first panel
h=prf.panels.L(prf.Nle-1)*2.5;

% guess second wake point 
guess= xLE+h*s;

xw=zeros(N,2); xw(1,:)=xLE;
sw=zeros(1,N);
lw=zeros(1,N-1);
%value of streamfunktion for guessing point

dy=30*abs(h*min(s(1),s(2)));

for i=1:N-1 % each node of wake
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
end

%tangential vektors

ew=transpose(xw(2:end,:)-xw(1:end-1,:)); ew=[ew(1,:)./lw; ew(2,:)./lw];
nw=[-ew(2,:);ew(1,:)];

% normal vector of node -> mean value of panel normal vectors
nwn=(nw(:,1:end-1)+nw(:,2:end))/2;
nwn=[nw(:,1),nwn, nw(:,end)];

LEstr.theta=atan2(ew(1,:),-ew(2,:));
LEstr.x=xw(:,1);
LEstr.y=xw(:,2);
LEstr.e=ew;
LEstr.n=nw;
LEstr.L=lw;
LEstr.s=sw;
LEstr.nn=nwn;
LEstr.N=N;


end

