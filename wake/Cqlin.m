function [ Cq ] = Cqlin(X1,X2,Y,r1,r2,lnr1,lnr2,t1,t2,cor,k1,k2, L )
%GRADPSIW calculates the source coefficients for velocity U=grad(Psi) \cdot n



% piecewise linear ansatz

% midpoint quantities
XM=(X1+X2)/2;
rM=XM.^2+Y.^2;
lnrM=log(rM)/2; 

sgn=sign(Y);
tM=atan2(sgn.*XM,sgn.*Y) + pi/2*(1-sgn);
tM=tM-cor;

n=length(XM(:,1));


% length of j-1 panel + j panel
Lm= [L(1), L(1:end-1) + L(2:end)]; 
hm=ones(n,1)*Lm;
% length of j panel
h=ones(n,1)*L;
% length of j+1 panel + j panel
Lp= [L(1:end-1) + L(2:end), L(end)]; 
hp=ones(n,1)*Lp;


 % Approximation of first distribution
pp= XM.*tM - X1.*t1 + Y.*(lnr1-lnrM);
pm= ( (X1+XM).*pp + r1.*t1 - rM.*tM + Y.*(XM-X1) )./(X1-XM);

dpm_dX1= t1 + (pp - pm)./(X1-XM);
dpm_dXM= tM + (pp + pm)./(X1-XM);
dpm_dY = ( (X1+XM).*(lnr1-lnrM) + 2*Y.*(t1-tM) )./(X1-XM) - 2;

gpp = k1.*(tM-t1) + k2.*(lnr1-lnrM);
gpm = k1.*(dpm_dX1 + dpm_dXM) + dpm_dY.*k2;

Cjm=  ( -gpp + gpm )./(4*pi*hm);
Cj = -(  gpp + gpm )./(4*pi*h);
Cjp=  (  gpp.*(1./h + 1./hm) + gpm.*(1./h - 1./hm) )/(4*pi);

Cj =[Cj ,zeros(n,1)];
Cjm=[Cjm(:,1) + Cjm(:,2),Cjm(:,3:end) ,zeros(n,2)];
Cjp=[zeros(n,1),Cjp];

C1= Cjm + Cj + Cjp;

% % Approximation of second distribution
pp= X2.*t2 - XM.*tM + Y.*(lnrM-lnr2);
pm= ( (X2+XM).*pp + rM.*tM - r2.*t2 + Y.*(X2-XM) )./(XM-X2);

dpm_dXM= tM + (pp - pm)./(XM-X2);
dpm_dX2= t2 + (pp - pm)./(XM-X2);
dpm_dY =  ((X2+XM).*(lnrM-lnr2) + 2*Y.*(tM-t2) )./(XM-X2) - 2;

gpp = k1.*(t2-tM) + k2.*(lnrM-lnr2);
gpm = k1.*(dpm_dXM + dpm_dX2) + dpm_dY.*k2;

Cj  = -(  gpp.*(1./h + 1./hp) + gpm.*(1./hp - 1./h) )/(4*pi);
Cjp =  (  gpp - gpm )./(4*pi*h);
Cjpp=  (  gpp + gpm )./(4*pi*hp);

Cj =[Cj ,zeros(n,1)];
Cjp=[zeros(n,1),Cjp];
Cjpp=[zeros(n,2), Cjpp(:,1:end-2) ,Cjpp(:,end-1)+Cjpp(:,end)];

C2= Cj + Cjp + Cjpp;


Cq=(C1+C2);


end

