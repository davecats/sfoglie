function [ Cg, Cq ] = GradPsiN( prf,wake )
%GRADPSIG  calculates the coefficients Cq_ij and Cg_ij for velocity U=grad(Psi) \cdot n for
%          i=1,..,Nw ;


ejx= prf.panels.e(1,:);
ejy= prf.panels.e(2,:);
nix= [wake.n(1,1); wake.n(1,:)'];
niy= [wake.n(2,1); wake.n(2,:)'];

% nix=wake.nn(1,:)';
% niy=wake.nn(2,:)';

% %k1=n1i*e1j + n2i*e2j =dX_dn
k1=  nix*ejx +  niy*ejy ;

% %k2=n2i*e1j - n1i*e2j
k2=  niy*ejx -  nix*ejy;


L=prf.panels.L;
NW=length(wake.x);

% airfoil part -> j=1,..,N

[X1,X2,Y]=GetLocalRelCoord(wake.x,wake.y,prf);


r1=X1.^2 + Y.^2; % r^2
r2=X2.^2 + Y.^2;

%Loadingpoint only for wake -> no singularity here
lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

% correct reflection on wake
sgn=sign(Y);
t1=atan2(sgn.*X1,sgn.*Y) + pi/2*(1-sgn);
t2=atan2(sgn.*X2,sgn.*Y) + pi/2*(1-sgn);


% circulation
%----------------------------------------------------------------

LM=ones(NW,1)*L;
pp= X1.*lnr1 - X2.*lnr2 - LM + Y.*(t1-t2);
pm= ( (X1+X2).*pp + r2.*lnr2 - r1.*lnr1 + 0.5*(X1.^2 - X2.^2) )./LM;

dpm_dX1 = ( (X2-X1).* lnr1 + pp - pm )./LM;
dpm_dX2 = ( (X2-X1).* lnr2 + pp + pm )./LM;
dpm_dY  = ( (X1+X2).* (t1-t2) + 2*Y.*(lnr2-lnr1) )./LM;

dpp_dn = k1.*(lnr1-lnr2) + k2.*(t1-t2);
dpm_dn = k1.*(dpm_dX1 + dpm_dX2) + k2.*dpm_dY;

C1= 1/(4*pi) * (dpp_dn - dpm_dn );
C2= 1/(4*pi) * (dpp_dn + dpm_dn );

% extract TE panel contribution
gTE=dpp_dn(:,end);

C1=[C1(:,1:end-1), zeros(NW,1)];
C2=[zeros(NW,1), C2(:,1:end-1)];

% Coefficient matrix for gammas
Cg=C1+C2;


% Trailing egde
qTE= k1(:,end).*( t2(:,end)-t1(:,end) ) + k2(:,end).*( lnr1(:,end) - lnr2(:,end) ) ;

TE= (prf.ScrossT*qTE - prf.SdotT*gTE)/(2*pi);

Cg(:,1)=Cg(:,1) + TE/2; 
Cg(:,end)=Cg(:,end) - TE/2; 


% source
%---------------------------

% also include wake part j=1,..,N +NW
[X1w,X2w,Yw]=WakeRelCoord(wake.x,wake.y,wake);



r1w=X1w.^2 + Yw.^2; % r^2
r2w=X2w.^2 + Yw.^2;


lnr1w=log(r1w)/2; 
lnr2w=log(r2w)/2;

sgn=sign(Yw);
t1w=atan2(sgn.*X1w,sgn.*Yw) + pi/2*(1-sgn);
t2w=atan2(sgn.*X2w,sgn.*Yw) + pi/2*(1-sgn);

% Filter singularity points out
S1= find(r1w<1e-11);
lnr1w(S1)=0;
t1w(S1)=0;
S2= find(r2w<1e-11);
lnr2w(S2)=0;
t2w(S2)=0;


% %k1=n1i*e1j + n2i*e2j =dX_dn
k1w=  nix*wake.e(1,:) +  niy*wake.e(2,:) ;

% %k2=n2i*e1j - n1i*e2j
k2w=  niy*wake.e(1,:) -  nix*wake.e(2,:);


% correct panel angle
cor= ones(NW,1)*prf.panels.theta;
t1=t1-cor;
t2=t2-cor;
corw= ones(NW,1)*wake.theta;
t1w=t1w-corw;
t2w=t2w-corw;

% piecewise linear ansatz
Cqwake= Cqlin(X1w,X2w,Yw,r1w,r2w,lnr1w,lnr2w,t1w,t2w,corw,k1w,k2w, wake.L );

% without TE -> TE source is part of gamma
CqFoil= Cqlin(X1(:,1:end-1),X2(:,1:end-1),Y(:,1:end-1),r1(:,1:end-1),r2(:,1:end-1),lnr1(:,1:end-1),lnr2(:,1:end-1),...
                t1(:,1:end-1),t2(:,1:end-1),cor(:,1:end-1),k1(:,1:end-1),k2(:,1:end-1), prf.panels.L(1:end-1) );

Cq=-[CqFoil, Cqwake];

% constant ansatz
%Cqt= ( k1.*(t2-t1) + k2.*(lnr1-lnr2) )/(2*pi);



end


