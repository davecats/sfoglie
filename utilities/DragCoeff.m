function cw = DragCoeff( T,H,Urel , mode)
%DRAGCOEFF calculates the Drag coefficient with the values on last wake node
%          mode=1: uses Formula of Eppler 2003
%          mode=2: uses Squire-Young formula 

if nargin==3
    mode=1;
end

kappa=0.5;
if mode ==1 % Epplers Formula    
    cw=2*Urel*T* (Urel-(1-Urel)*(1-kappa)*H)/(1+kappa*(1-Urel));
else % Squire-Young-formula
    cw=2*T*Urel^( 3 + kappa*(H-1) );
end


end

