function [ Cq ] = Cqlin(X1,X2,Y,r1,r2,lnr1,lnr2,t1,t2,cor,k1,k2, L )
%Cqlin calculates the source coefficients for velocity U=grad(Psi) \cdot n
%in case of a linear ansatz between the panels midpoints. 
% Input:
% X1_ij= ( xi_i - x1_j ) \cdot e_j                   distance of Loading point and start node in x-direction of the local coord system
% X2_ij= ( xi_i - x2_j ) \cdot e_j = L_j - X1_ij
% Y_ij = ( xi - x1 ) \cdot n                         distance of Loading point to panel in y-direction of the local coord system
% r1_ij= ||xi_i - x1_j||^2 = X1^2 + Y^2              square of the norm of relativ vector xi - x1
% t1_ij= atan2(X1,Y)                                 Angle of relativ vector to local y-axis
% k1_ij=n1_i*e1_j + n2_i*e2_j = dX_i /dn_j          derivate of x-coord of i-th loading point rel vec in direction of j-th panel normal vec   
% k2_ij=n2_i*e1_j - n1_i*e2_j = dY_i /dn_j          derivate of y-coord of i-th loading point rel vec in direction of j-th panel normal vec

% piecewise linear ansatz

% midpoint quantities
XM=(X1+X2)/2;
rM=XM.^2+Y.^2;
lnrM=log(rM)/2; 

% adjust atan branching if necessary
sgn=sign(Y);
tM=atan2(sgn.*XM,sgn.*Y) + pi/2*(1-sgn);
% add the panel angle
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


 % Approximation of the distribution for the first panel half
pp= XM.*tM - X1.*t1 + Y.*(lnr1-lnrM);
pm= ( (X1+XM).*pp + r1.*t1 - rM.*tM + Y.*(XM-X1) )./(X1-XM);

% Derivates for the gradient
dpm_dX1= t1 + (pp - pm)./(X1-XM);
dpm_dXM= tM + (pp + pm)./(X1-XM);
dpm_dY = ( (X1+XM).*(lnr1-lnrM) + 2*Y.*(t1-tM) )./(X1-XM) - 2;

% inner product of gradient and normal vector grad(Psi) \cdot n
gpp = k1.*(tM-t1) + k2.*(lnr1-lnrM);
gpm = k1.*(dpm_dX1 + dpm_dXM) + dpm_dY.*k2;

Cjm=  ( -gpp + gpm )./(4*pi*hm);
Cj = -(  gpp + gpm )./(4*pi*h);
Cjp=  (  gpp.*(1./h + 1./hm) + gpm.*(1./h - 1./hm) )/(4*pi);

Cj =[Cj ,zeros(n,1)];
Cjm=[Cjm(:,1) + Cjm(:,2),Cjm(:,3:end) ,zeros(n,2)];
Cjp=[zeros(n,1),Cjp];

C1= Cjm + Cj + Cjp;

% % Approximation of distribution for the second panel half
pp= X2.*t2 - XM.*tM + Y.*(lnrM-lnr2);
pm= ( (X2+XM).*pp + rM.*tM - r2.*t2 + Y.*(X2-XM) )./(XM-X2);

% Derivates for the gradient
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

% Add coefficients of both panel half for the total matrix
Cq=(C1+C2);


end

