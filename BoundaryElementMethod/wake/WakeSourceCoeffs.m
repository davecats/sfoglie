function [ Bwake, Cqwake ] = WakeSourceCoeffs( wake,prf )
%WAKESOURCECOEFFS   calculates the influence of wake sources to all equations:
%                      Bwake=Bij,  i=1,..,N j=N+1,..,N+NW
%                           -> influence of wake sources on airfoil node Eqs
%                      Cqwake=Cij,  i=N+1,..,N+NW j=N+1,..,N+NW
%                           -> influence of wake sources on wake node Eqs


NW=length(wake.x);
N=prf.N;
% airfoil velocity
%---------------------------

[X1,X2,Y]=WakeRelCoord(transpose(prf.nodes.X),transpose(prf.nodes.Y),wake);

r1=X1.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

t1=atan2(X1,Y);
t2=atan2(X2,Y);

% no singularities to filter

% constant ansatz -> only q_j appears in integral for j-th panel
psi = (1/(2*pi))*( -X1.*t1 +X2.*t2 + Y.*(lnr1-lnr2) );
Bwake=[psi, zeros(N,1)]; % no source contribution on last wake node


% Wake velocity
%---------------------------

[X1,X2,Y]=WakeRelCoord(wake.x,wake.y,wake);

r1=X1.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

Y=-Y;
t1=atan2(X1,Y);
t2=atan2(X2,Y);

% Filter singularity points out
S1= find(r1<1e-11);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-11);
lnr2(S2)=0;
t2(S2)=0;

% k1=e1*n1+e2*n2
k1= wake.e(1,:).*wake.n(1,:) +  wake.e(2,:).*wake.n(2,:);
k1=ones(NW,1)*k1;
% k2=e1*n2-e2*n1 = - (n1^2 + n2^2)
k2= wake.e(1,:).*wake.n(2,:) -  wake.e(2,:).*wake.n(1,:);
k2=ones(NW,1)*k2;


Cqwake=   k1.*(t2-t1) + k2.*(lnr1-lnr2)/(4*pi);
Cqwake=[Cqwake, zeros(NW,1)];%last value doesn´t appear for linear ansatz


end

