function [wake] = GetWakeStreamline( fld,prf,NW ,grading)
%GETWAKESTREAMLINE  calculates the noedes on the wake streamline. 
%                   Finds points where psi-psi0 goes to zero -> Sekantenmethode
%                   NW: number of wake nodes
%                   grading: default set to 1.18


if nargin==3
grading=1.18;
end

%        calculate TE midpoint as first wake node
x1=[prf.panels.X(end,2) prf.panels.Y(end,2)]; 
xN=[prf.panels.X(end,1) prf.panels.Y(end,1)];

xTE= (x1+xN)/2; % Midpoint of TE

%psi0=fld.psi0;
psi0=evaluateInviscFieldSol(xTE,fld,prf,'without TE');


%        calculate first wake tangent vector
%           -> mean of tangent vectors of first and N-th panel 
e1=prf.panels.e(1,:);
eN=prf.panels.e(end-1,:);
s=0.5*(eN-e1);
%if s(1)<0; s=-s; end;




% starting step size of wake equals the Panel length of first panel
h=prf.panels.L(1);

% guess second wake point
guess= xTE+2.5*h*s;

xw=zeros(NW,2); xw(1,:)=xTE;
sw=zeros(NW,1);
lw=zeros(NW-1,1);
%value of streamfunktion for guessing point

dy=30*abs(h*s(2));

for i=1:NW-1 % each node of wake
    psig=evaluateInviscFieldSol(guess,fld,prf,'without TE');
    found=false;k=1;
    
    % find point on wake streamline -> psi=psi0
    
    % Einschlussintervall finden
    %------------------------------------------------
    while found==false && k<20 % prevent endless loop
        psip=evaluateInviscFieldSol(guess-[0, k*dy],fld,prf,'without TE');
        d1=psig-psi0;
        d2=psip-psi0;
        if sign(d1)*d2 < 0 % Einschlussintervall gefunden
          psi1=psig;
          psi2=psip;
          found=true;
          y2=guess(2)-k*dy;
        else
          psim=evaluateInviscFieldSol(guess+[0, k*dy],fld,prf,'without TE');  
          d2=psim-psi0;  
          if sign(d1)*d2 < 0 % Einschlussintervall gefunden
            psi1=psig;
            psi2=psim;  
            found=true;
            y2=guess(2)+k*dy;
          end
        end
        k=k+1;   
    end
    %------------------------------------------------
    y1=guess(2);
    dn=1; l=0;
    % Sekanten Methode
    %------------------------------------------------
    res=1e-11; % residuum
    while  abs(dn)>res && l<40 % prevent endless loop
        d1=psi1-psi0;
        d2=psi2-psi0;
        yn=(abs(d1)*y2+abs(d2)*y1) /(abs(d1)+abs(d2));
        new=[guess(1), yn];
        psin=evaluateInviscFieldSol(new,fld,prf,'without TE');
        dn=psin-psi0;
        if sign(dn)*d1<0
           psi2=psin;
           y2=yn;
        else
           psi1=psin;
           y1=yn; 
        end
        l=l+1;
    end
    %------------------------------------------------
    xw(i+1,:)=[guess(1), yn];
    guess=xw(i+1,:)+grading*(xw(i+1,:)-xw(i,:));
    lw(i)=norm(xw(i+1,:)-xw(i,:));
    sw(i+1)=sw(i)+lw(i);
    
end

%tangential vektors

ew=(xw(2:end,:)-xw(1:end-1,:));ew=[ew(:,1)./lw, ew(:,2)./lw];
nw=[ew(:,2),-ew(:,1)];

% normal vector of node -> mean value of panel normal vectors
nwn=(nw(1:end-1,:)+nw(2:end,:))/2;nwn=[nw(1,:);nwn; nw(end,:)];


wake.x=xw(:,1);
wake.y=xw(:,2);
wake.e=ew;
wake.n=nw;
wake.L=lw;
wake.s=sw;
wake.nn=nwn;
%{
for i=1:NW-1
    st= false;
    psinew=evaluateInviscFieldSol(guess,fld,prf); %streamfunktion value at guessed point
    while st==false
        dy=14*abs(h*s(2));
        psi2=evaluateInviscFieldSol(guess+[0, dy],fld,prf);
        psi3=evaluateInviscFieldSol(guess-[0, dy],fld,prf);
        d1=psinew-psi0;
        d2=psi2-psi0;
        d3=psi3-psi0;
        if abs(d2)<abs(d1) %get correction direction
            guess=guess+[0, dy];
            psinew=psi2;
        elseif abs(d3)<abs(d1)
            guess=guess-[0, dy];
            psinew=psi3;
        else
            st=true;
        end
    end
    xw(i+1,:)=guess;
    guess=xw(i+1,:)+grading*(xw(i+1,:)-xw(i,:));
    sw(i+1)=sw(i)+norm(xw(i+1,:)-xw(i,:));
end
%}

end

