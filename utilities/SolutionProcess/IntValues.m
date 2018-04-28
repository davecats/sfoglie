function [CL,Cdrag,Cnu] = IntValues( Cp, tau, s, e, n,alfa,simpson,Vb)
%INTVALUES numerically integrates the values CL and Cdrag over the profile
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


% add blowing contribution
if nargin==8
   Cp = Cp + sign(Vb(1:length(Cp))).*Vb(1:length(Cp)).^2;   
end

% pressure forces of each panel in x and y direction 
fp_x= Cp.*n(:,1);
fp_y= Cp.*n(:,2);


% integrate with simpson law
FR_x= 2*NumInt(fR_x,s,mode); % factor 2 because of tau=cf/2
FR_y= 2*NumInt(fR_y,s,mode);
FP_x= NumInt(fp_x,s,mode);
FP_y= NumInt(fp_y,s,mode);


% total force
Fx=FR_x + FP_x;
Fy=FR_y + FP_y;


% Transform to coordinate system in stream direction
FX=  Fx*cos(alfa) + Fy*sin(alfa);
FY= -Fx*sin(alfa) + Fy*cos(alfa);

CL=FY;
Cdrag=abs(FX);
Cnu= abs( FR_x*cos(alfa) );






end


% figure
% hold on
% plot(s,fR_x)
% plot(s,fp_x)
% legend('Reibungsanteil','Druckanteil')
% 






