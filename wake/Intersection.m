function [ new, psin ] = Intersection( fld,prf,guess,h,psi0 )
%INTERSECTION Intersectionmethod to find wake streamline -> in case Newton
%             method does not converge



itmax=30;
psig=evaluateInviscFieldSol(guess,fld,prf);
found=false;k=1;

% find point on wake streamline -> psi=psi0

% find intersection intervall
%------------------------------------------------

dy=0.1*h;


while found==false && k<itmax % prevent endless loop
    psip = evaluateInviscFieldSol(guess+[0, k*dy],fld,prf);
    d1=psig-psi0;
    d2=psip-psi0;
    if sign(d1)*d2 < 0 % intersection intervall found
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

if ~found
   error('No intersectionintervall found')
   return 
end
%------------------------------------------------
y1=guess(2);
dn=1; l=0;
% Secant method
%------------------------------------------------
res=1e-10; % residuum
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



end

