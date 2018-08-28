function [CM,MR,Mp ] = getCM(prf,sol,flo,refPoint,simpson)
%GETCM  calculates the momentum coefficient for the profile 
%       refPoint: reference point for the momentum
%       simpson=true -> numerical integration with simpson law
%       simpson=false-> numerical integration with mitpoint rule

if nargin<4 % ref point not committed
   refPoint=[0.25,0]*prf.c;
   disp('Momentum coefficient reference point unspecified -> set to x=0.25c y=0')
end
if nargin<5 || simpson
    mode=1;
else
    mode=2;
end

% get the stress vectors 
[~,fR,fp] = StressVector(sol.Cp(1:prf.N),2*sol.tau(1:prf.N),prf.nodes.n',sol.Vb/flo.Uinfty );

% relative vector of each node to reference point
relX=prf.nodes.X'-refPoint(1);
relY=prf.nodes.Y'-refPoint(2);

% momentum density   
% cross product m=  f2 relX - f1 relY
mR=fR(:,2).*relX - fR(:,1).*relY;
mp=fp(:,2).*relX - fp(:,1).*relY;

% numerical integration
MR=NumInt(mR(:,1),prf.s,mode);
Mp=NumInt(mp(:,1),prf.s,mode);

% add pressure and friction part
CM=MR+Mp;



end


