function [CL,Cw,CwR,Cwp] = IntValues( Cp, tau, s, n,alfa,simpson,Vb)
%INTVALUES numerically integrates the shear and pressure forces on the
%          profile to gain lift and drag coefficient
%          input:
%          Cp  : pressure coefficient
%          tau : wall shear stress (normalized with rho Uinfty^2 -> without factor 2 !!!) 
%          s   : profile arc length vector
%          n   : normal vector of nodes
%          alfa: angle of attack 
%          Vb  : wall normal velocity normalized with Uinfty! (default is zero) 
%          simpson= true -> uses simpson law for integration (default)
%          simpson= false-> uses mid point rule 
%           
%          output:
%          CL : lift coefficient
%          Cw : drag coefficient
%          CwR: friction drag coefficient
%          Cwp: pressure drag coefficient 


if nargin==5 || simpson
    mode=1;
else
    mode=2;
end

% calculate stresses
[~,fR,fp]= StressVector(Cp,tau,n,Vb );


% numerical integration
FR_x= 2*NumInt(fR(:,1),s,mode); % factor 2 because of tau=cf/2
FR_y= 2*NumInt(fR(:,2),s,mode);
FP_x= NumInt(fp(:,1),s,mode);
FP_y= NumInt(fp(:,2),s,mode);


% total force (add upp shear and pressure contribution)
Fx=FR_x + FP_x;
Fy=FR_y + FP_y;


% Transform to coordinate system in stream direction
FX=  Fx*cos(alfa) + Fy*sin(alfa);
FY= -Fx*sin(alfa) + Fy*cos(alfa);

% output values
CL=FY;
Cw=FX;
CwR=  FR_x*cos(alfa);
Cwp=Cw-CwR;


end







