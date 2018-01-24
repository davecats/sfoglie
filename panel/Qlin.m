function B = Qlin( prf, wake)
%QLIN  calculates source coefficients B with a two piece linear ansatz
%      if only prf is comitted -> Bij, i=1,..,N ; ij=1,..,N
%      if wake is comitted     -> Bij, i=1,..,N ; ij=N+1,..,N+NW


if nargin==1 % Contribution of airfoil nodes
L=prf.panels.L(1:end-1); % panel length without TE panel
N=length(L); % number of panels (without TE panel), N+1 number of nodes


% panel angle
panAng=prf.panels.theta(1:end-1);

% node coordinates 
[X1,X2,Y]=GetLocalRelCoord(transpose(prf.nodes.X),transpose(prf.nodes.Y),prf,'noTE');
else % Contribution of wake nodes
L=wake.L;
N=length(prf.nodes.X)-1;
panAng=wake.theta;

[X1,X2,Y]=WakeRelCoord(transpose(prf.nodes.X),transpose(prf.nodes.Y),wake);
end

% panel midpoint
XM=(X1+X2)/2;


r1=X1.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;
rM=XM.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;
lnrM=log(rM)/2; 

t1=atan2(X1,Y); 
t2=atan2(X2,Y); 
tM=atan2(XM,Y);



% t1=atan(X1./Y); 
% t2=atan(X2./Y); 
% tM=atan(XM./Y);

% sgn=ones(N+1,N+1);
% sgn(1,2:end)=sign(Y(1,2:end));
% sgn(N+1,1:end-1)=sign(Y(N+1,1:end-1));
% t1=atan2(sgn.*X1,sgn.*Y) + (ones(N+1,N+1) - sgn)*pi/2;
% t2=atan2(sgn.*X2,sgn.*Y) + (ones(N+1,N+1) - sgn)*pi/2;
% tM=atan2(sgn.*XM,sgn.*Y) + (ones(N+1,N+1) - sgn)*pi/2;

% Filter singularity points out
S1= find(r1<1e-11);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-11);
lnr2(S2)=0;
t2(S2)=0;


% a=load('./XFoilWerte/t1.txt');
% diff=a(:,1:end-1)-t1;
% max(max(abs(a(:,1:end-1)-t1)))

cor=ones(N+1,1)*panAng;
t1=t1-cor; 
t2=t2-cor; 
tM=tM-cor;

% length of j-1 panel + j panel
Lm= [L(1), L(1:end-1) + L(2:end)]; 
hm=ones(N+1,1)*Lm;
% length of j panel
h=ones(N+1,1)*L;
% length of j+1 panel + j panel
Lp= [L(1:end-1)+L(2:end), L(end)]; 
hp=ones(N+1,1)*Lp;


% % Approximation of first distribution
pp= XM.*tM - X1.*t1 + Y.*(lnr1-lnrM);
pm= ( (X1+XM).*pp +r1.*t1 - rM.*tM + Y.*(XM-X1) )./(X1-XM);

% pDX1= (-(X1+XM).*t1 + pp + 2*X1.*t1 - pm)./(X1-XM);
% pDXM= ( (X1+XM).*tM + pp - 2*XM.*tM + pm)./(X1-XM);
% pDY = ( (X1+XM).*(lnr1-lnrM) + 2*(XM-X1 + Y.*(t1-tM)) )./(X1-XM);

Cjm=  ( -pp + pm )./(4*pi*hm);
Cj = -(  pp + pm )./(4*pi*h);
Cjp=  (  pp.*(1./h+1./hm) + pm.*(1./h-1./hm) )/(4*pi);


Cj =[Cj ,zeros(N+1,1)];
Cjm=[Cjm(:,1)+Cjm(:,2),Cjm(:,3:end) ,zeros(N+1,2)];
Cjp=[zeros(N+1,1),Cjp];

B1= Cjm + Cj + Cjp;

% % Approximation of second distribution
pp= X2.*t2 - XM.*tM + Y.*(lnrM-lnr2);
pm= ( (X2+XM).*pp +rM.*tM - r2.*t2 + Y.*(X2-XM) )./(XM-X2);

% pDX1= (-(X2+XM).*tM + pp + 2*XM.*tM - pm)./(XM-X2);
% pDXM= ( (X2+XM).*t2 + pp - 2*X2.*t2 + pm)./(XM-X2);
% pDY =  ((X2+XM).*(lnrM-lnr2) + 2*(XM-X1 + Y.*(t1-tM)) )./(XM-X2);

Cj  = -(  pp.*(1./h+1./hp) + pm.*(1./hp-1./h) )/(4*pi);
Cjp =  (  pp - pm )./(4*pi*h);
Cjpp=  (  pp + pm )./(4*pi*hp);

Cj =[Cj ,zeros(N+1,1)];
Cjp=[zeros(N+1,1),Cjp];
Cjpp=[zeros(N+1,2), Cjpp(:,1:end-2) ,Cjpp(:,end-1)+Cjpp(:,end)];

B2= Cj + Cjp + Cjpp;


B=(B1+B2);



end

