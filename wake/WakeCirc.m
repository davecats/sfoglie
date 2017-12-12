function [ c1, c2 ] = WakeCirc( extrema,point,ni )
%WAKECIRC calculates the influence of the circulation du_dgamma for the wake velocity


t=extrema(2,:)-extrema(1,:); % tangent vector of the panel
L=norm(t); t=t/L; 
n=[t(2) -t(1)]; %t=t; % normal vector of the panel
X1=sum((point-extrema(1,:)).*t);  %Delta x from first node to the loading point (x,y)
X2=sum((point-extrema(2,:)).*t); %Delta x from second node to the loading point (x,y)
Y=sum((point-extrema(1,:)).*n);  %Delta y of the panel to the loading point (x,y)
r1=norm(point-extrema(1,:));     %distance from first node to the loading point (x,y)
r2=norm(point-extrema(2,:));     %distance from second node to the loading point (x,y)


 sgn=sign(Y);
 % Filtering the singularity on bouandary points 
%   -> if (x,y) lays in one of the nodes of the panel
 if r1<1e-12; lnr1=0;  t1=0; else lnr1=log(r1); t1=atan2(sgn*X1,sgn*Y)+pi*(1-sgn)/2;  end
 if r2<1e-12; lnr2=0;  t2=0; else lnr2=log(r2); t2=atan2(sgn*X2,sgn*Y)+pi*(1-sgn)/2; end

 
 pp= X1*lnr1 - X2*lnr2 + X2 - X1 + Y*(t1-t2);
 pm= ((X1+X2)*pp + (r2*lnr2-r1*lnr1 + X1*X1-X2*X2))/(X1-X2);
 
 u1= ((X2-X1)*lnr1  +pp -pm)/(X1-X2);
 u2= ((X2-X1)*lnr2  +pp +pm)/(X1-X2);
 u3= ((X2-X1)*(t1-t2) -Y*(lnr1-lnr2))/(X1-X2);

 xN= sum(t.*ni);
 yN= t(1)*ni(2)-t(2)*ni(1);   

 up=(lnr1-lnr2)*xN+yN*(t1-t2);
 um=(u1-u2)*xN+yN*u3;
 
 c1=(up-um)/(4*pi);
 c2=(up+um)/(4*pi);
 
end

