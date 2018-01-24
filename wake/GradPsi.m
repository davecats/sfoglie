function [ Cg, Cq ] = GradPsi( prf,wake )
%GRADPSIG  calculates the coefficients for velocity U=grad(Psi) \cdot n for
%          i=1,..,Nw ; j=1,..,N -> airfoil influence on wake node equations
%           -> U_i= Cg_ij gamma_j + Cq_ij q_j


n=prf.panels.n(1:end,:); 
e=prf.panels.e(1:end,:);
L=prf.panels.L(1:end,:);


NW=length(wake.x);

[X1,X2,Y]=GetLocalRelCoord(wake.x,wake.y,prf);


r1=X1.^2+Y.^2; % r^2
r2=X2.^2+Y.^2;

%Loadingpoint only for wake -> no singularity here
lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

t1=atan2(X1,Y);
t2=atan2(X2,Y);
 
% circulation
%---------------------------
% k1=e1*n1+e2*n2 =dX_dn
k1= e(1,:).*n(1,:) +  e(2,:).*n(2,:);
k1=ones(NW,1)*k1;
% k2=e1*n2-e2*n1 =  (n1^2 + n2^2)
k2= e(1,:).*n(2,:) -  e(2,:).*n(1,:) ;
k2=ones(NW,1)*k2;



LM=ones(NW,1)*L;
A=( (k2.*Y-k1.*X2).*(lnr2-lnr1) - (k1.*Y+k2.*X2).*(t2-t1) - k1.*LM   )./(2*pi*LM);
GTE=A(:,end);% extract trailing edge panel
A=[A(:,1:end-1), zeros(NW,1)]; % coefficient for j-th node in j-th panel -> last node doesn´t appear  

B=( (k2.*Y-k1.*X1).*(lnr1-lnr2) + (k1.*Y+k2.*X1).*(t2-t1) + k1.*LM   )./(2*pi*LM);
GTE=(GTE+B(:,end));% extract trailing edge panel
B=[zeros(NW,1),B(:,1:end-1)]; % coefficient for j+1-th node in j-th panel -> first node doesn´t appear 

Cg=A+B;


% source
%---------------------------

Cq= ( k1.*(t2-t1) + k2.*(lnr1-lnr2) )/(2*pi);%  - X1.*(k1.*Y-k2.*X1)./r1 + X2.*(k1.*Y-k2.*X2)./r2;

STE=Cq(:,end);
Cq= [Cq(:,1:end-1), zeros(NW,1) ]; %last value doesn´t appear for linear ansatz

% Add Trailing edge influence to matrix A
TE= prf.SdotT*GTE + prf.ScrossT*STE;
%Cg(:,1)=Cg(:,1)+TE; % with kutta
 Cg(:,1)=Cg(:,1)+TE/2; 
 Cg(:,end)=Cg(:,end)-TE/2;    
end

