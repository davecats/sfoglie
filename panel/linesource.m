function [psi] = linesource(extrema,point)
%LINESOURCE calculates the coefficients of the source q_j for the
%           j-th panel (constant ansatz functions for solution).          
%           -> each node appears only in the integral of one panel.
%
%           input: extrema - Matrix with panel startpoint and end point[X_1 , Y_1, X_2, Y_2]
%                  point -   loading point = the  i-th node    


t=extrema(2,:)-extrema(1,:); % tangent vector of the panel
L=norm(t); t=t/L; 
n=[t(2) -t(1)]; %t=t; % normal vector of the panel
X=sum((point-extrema(1,:)).*t);  %Delta x from first node to the loading point (x,y)
X2=sum((point-extrema(2,:)).*t); %Delta x from second node to the loading point (x,y)
Y=sum((point-extrema(1,:)).*n);  %Delta y of the panel to the loading point (x,y)
r1=norm(point-extrema(1,:));     %distance from first node to the loading point (x,y)
r2=norm(point-extrema(2,:));     %distance from second node to the loading point (x,y)


 sgn=sign(X);
 sgn2=sign(X2);
 % Filtering the singularity on bouandary points 
 %   -> if (x,y) lays in one of the nodes of the panel
 if r1<1e-10; lnr1=0; t1=0;  else lnr1=log(r1);  t1=atan2(sgn*X,sgn*Y)+pi/2*(1-sgn); end
 if r2<1e-10; lnr2=0; t2=0;  else lnr2=log(r2);  t2=atan2(sgn2*X2,sgn2*Y)+pi/2*(1-sgn2); end
t1=t1- InvAngle(atan2(-t(1),t(2) )+pi);
t2=t2- InvAngle(atan2(-t(1),t(2) )+pi);

 
psi = (1/(2*pi))*( -X*t1 +X2*t2 + Y*(lnr1-lnr2) );

% Alternative formulation 
%
% if r1<1e-10; lnr1=0; t1=0; else t1=atan2(Y,X);  lnr1=log(r1); end
% if r2<1e-10; lnr2=0; t2=0; else t2=atan2(Y,X2); lnr2=log(r2); end
% 
% 
% if t1<0; t1=2*pi+t1; end; t1=t1-pi/2; 
% if t1<0; t1=3/2*pi+(pi/2+t1); end
% 
% if t2<0; t2=2*pi+t2; end; t2=t2-pi/2;   
% if t2<0; t2=3/2*pi+(pi/2+t2); end
% 
% psi = (1/(2*pi))*( X*t1 - X2*t2 + Y*(lnr1-lnr2) );


