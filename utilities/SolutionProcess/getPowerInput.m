function [Pb] = getPowerInput(Cp,vw,Uinfty,s, IntMode )
%GETPOWERINPUT calculates the normalized power input required for uniform
%blowing P/(0.5*rho*Uinfty^3): 
%           int Cp*vw/Uinfty + (vw/Uinfty)^3 ds
% Cp: pressure coefficient
% vw: blowing velocity
% Uinfty: inflow velocity
% s: arclengthvector of the airfoil
% IntMode: mode for numerical integration (1 Simpson law, 2 mid point law)

if nargin==4; IntMode=1; end % default -> uses Simpson law for integration

n=length(s);

%dimensionless
vwD= vw(1:n)/Uinfty;

% each node
Pbdist = Cp(1:n).*vwD + vwD.^3;

%figure
%plot(s,Pbdist)
% prevent a power gain on suction side
Pbdist(Pbdist<0)=0;


% total power input -> integration over the whole distribution
Pb = NumInt(Pbdist,s,IntMode);

end

