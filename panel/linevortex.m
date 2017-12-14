function [psi1, psi2] = linevortex(extrema,point)
%LINEVORTEX calculates the coefficients of gamma_j and gamma_j+1 for the
%           j-th panel (linear ansatz functions for solution).          
%           -> each node appears in the integrals of two panels.
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


 sgn=sign(Y);
 % Filtering the singularity on bouandary points 
%   -> if (x,y) lays in one of the nodes of the panel
 if r1<1e-12; lnr1=0;  t1=0; else lnr1=log(r1); t1=atan2(sgn*X,sgn*Y)+pi*(1-sgn)/2;  end
 if r2<1e-12; lnr2=0;  t2=0; else lnr2=log(r2); t2=atan2(sgn*X2,sgn*Y)+pi*(1-sgn)/2; end
% Neue Koeffizienten ---------------------------------------
%coefficient for the first and second node of the panel
psi1=-(1/(4*pi*L))*( (X*(X-2*L)-Y^2)*lnr1 - ((L-X)^2-Y^2)*lnr2 ...
                     + 2*Y*(L-X)*(t2-t1) +L*(3*L/2-X) );
psi2=-(1/(4*pi*L))*( -1*(X^2-Y^2)*lnr1 +(X^2-Y^2-L^2)*lnr2 ...
                     + 2*Y*X*(t2-t1) +(1/2)*L*(L+2*X));



% if r1<1e-12; lnr1=0;  t1=0; else lnr1=log(r1); t1=atan2(Y,X);  end
% if r2<1e-12; lnr2=0;  t2=0; else lnr2=log(r2); t2=atan2(Y,X2); end
%if Y<1e-9; t1=0; t2=pi; end
% if abs(t1-t2)>pi+1e-5
%     t1*180/pi
%     t2*180/pi
%     error('eeeee'); 
% end

% psi1 = (1/(2*pi*L))*( -0.5*(X*(X-2*L)-Y^2)*lnr1  ...
%                       +0.5*((L-X)^2-Y^2)*lnr2 ...
%                       +(L-X)*Y*(t2-t1) ...
%                       -0.5*(L-X)*L );
% psi2 = (1/(2*pi*L))*( +0.5*(X^2-Y^2)*lnr1  ...
%                       -0.5*(X^2-Y^2-L^2)*lnr2 ...
%                       +X*Y*(t2-t1) ...
%                       -0.5*X*L );





                 
                  