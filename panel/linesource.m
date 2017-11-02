function [psi] = linesource(extrema,point)

t=extrema(2,:)-extrema(1,:); 
L=norm(t); t=t/L; 
n=[t(2) -t(1)]; t=t; 
X=sum((point-extrema(1,:)).*t);
X2=sum((point-extrema(2,:)).*t); 
Y=sum((point-extrema(1,:)).*n); 
r1=norm(point-extrema(1,:)); 
r2=norm(point-extrema(2,:)); 
if r1<1e-10; lnr1=0; t1=0; else; t1=atan2(Y,X);  lnr1=log(r1); end
if r2<1e-10; lnr2=0; t2=0; else; t2=atan2(Y,X2); lnr2=log(r2); end
% 
% if t1<0; t1=2*pi+t1; end; t1=t1-pi;
% if t2<0; t2=2*pi+t2; end; t2=t2-pi; 

%if t1<0; t1=3/2*pi+(pi/2+t1); end

if t1<0; t1=2*pi+t1; end; t1=t1-pi/2; 
if t1<0; t1=3/2*pi+(pi/2+t1); end

if t2<0; t2=2*pi+t2; end; t2=t2-pi/2;   
if t2<0; t2=3/2*pi+(pi/2+t2); end

psi = (1/(2*pi))*( X*t1 - X2*t2 + Y*(lnr1-lnr2) );

% 
% 
% z1 = (X+1i*Y);
% z2 = (X-L+1i*Y);
% if abs(z1)<1e-9; lnz1=0; else; lnz1=log(z1); end
% if abs(z2)<1e-9; lnz2=0; else; lnz2=log(z2); end
% F = -(z2*lnz2-z2-z1*lnz1+z1)/(2*pi*L); 
% psi=imag(F); 
