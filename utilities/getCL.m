function [CL ] = getCL(prf,U,mode)
%GETCL  calculates the lift coefficient for the profile with given boundary
%       edge tangential velocity 
%       mode=1: midpoint integration rule ( second order accurate) ->  default
%       mode=2: Simpson integration rule ( fourth order accurate)

if nargin==2
    mode=1;
end

CP= 1- (U(1:prf.N)/prf.Uinfty).^2;

% add TE CP
CP=[CP;CP(1)]; 

% panel midpoint value
CPM= 0.5* (CP(2:end)+CP(1:end-1));

% force distributions of panel midpoints
fpx = CPM.*prf.panels.n(1,:)';
fpy = CPM.*prf.panels.n(2,:)';

if mode==1 % midpointrule
   Fpx= sum( fpx.*prf.panels.L' );
   Fpy= sum( fpy.*prf.panels.L' );
else % simpson rule
   s=[prf.s ,prf.s(end)+prf.panels.L(end)];
   s= 0.5*(s(2:end)+s(1:end-1));
   
   Fpx=NumInt(fpx,s);
   Fpy=NumInt(fpy,s);   
end

% transformation in streamwise coordinate System
%FX=  Fpx*cos(prf.alfa) + Fpy*sin(prf.alfa);
FY= -Fpx*sin(prf.alfa) + Fpy*cos(prf.alfa);

CL=FY;


end


