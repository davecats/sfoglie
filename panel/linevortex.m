function [psi1, psi2] = linevortex(extrema,point)

t=extrema(2,:)-extrema(1,:); 
L=norm(t); t=t/L; 
n=[t(2) -t(1)]; t=t; 
X=sum((point-extrema(1,:)).*t);
X2=sum((point-extrema(2,:)).*t); 
Y=sum((point-extrema(1,:)).*n); 
r1=norm(point-extrema(1,:)); 
r2=norm(point-extrema(2,:)); 
if r1<1e-12; lnr1=0;  t1=0; else; lnr1=log(r1); t1=atan2(Y,X);  end
if r2<1e-12; lnr2=0;  t2=0; else; lnr2=log(r2); t2=atan2(Y,X2); end
%if Y<1e-9; t1=0; t2=pi; end
% if abs(t1-t2)>pi+1e-5
%     t1*180/pi
%     t2*180/pi
%     error('eeeee'); 
% end


psi1 = (1/(2*pi*L))*( -0.5*(X*(X-2*L)-Y^2)*lnr1  ...
                      +0.5*((L-X)^2-Y^2)*lnr2 ...
                      +(L-X)*Y*(t2-t1) ...
                      -0.5*(L-X)*L );
psi2 = (1/(2*pi*L))*( +0.5*(X^2-Y^2)*lnr1  ...
                      -0.5*(X^2-Y^2-L^2)*lnr2 ...
                      +X*Y*(t2-t1) ...
                      -0.5*X*L );
                 
                  