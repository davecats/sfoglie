function [CL,Cdrag,Cnu] = IntValues( Cp, tau, s, e, n,alfa,simpson,Vb)
%INTVALUES numerically integrates thevalues CL and Cdrag over the Profil
%           simpson=true: uses simpson law for integration (default)
%           simpson=false: uses midpoint rule
%           Vb-> adds contribution of blowing velocity

if nargin==6 || simpson
    mode=1;
else
    mode=2;
end

Nle=find(e(:,1)>0,1);

% shear forces of each panel in x and y direction 
fR_x= tau.*e(:,1);
fR_y= tau.*e(:,2);

% correct sign for suction side -> force is in opposit direction of panel tangent vector
fR_x(1:Nle-1)=-fR_x(1:Nle-1);
fR_y(1:Nle-1)=-fR_y(1:Nle-1);

if nargin==8
   Cp = Cp + 0.5*Vb(1:length(Cp)).^2;   
end

% pressure forces of each panel in x and y direction 
fp_x= Cp.*n(:,1);
fp_y= Cp.*n(:,2);

% add blowing contribution

    
% integrate with simpson law
FR_x= -NumInt(fR_x,s,mode);
FR_y= -NumInt(fR_y,s,mode);
FP_x= NumInt(fp_x,s,mode);
FP_y= NumInt(fp_y,s,mode);

% total force
Fx=2*FR_x + FP_x;
Fy=2*FR_y + FP_y;




% Transform to coordinate system in stream direction
FX=  Fx*cos(alfa) + Fy*sin(alfa);
FY= -Fx*sin(alfa) + Fy*cos(alfa);

CL=FY;
Cdrag=abs(FX);
Cnu= 2*sqrt(FR_x^2 + FR_y^2);






end

