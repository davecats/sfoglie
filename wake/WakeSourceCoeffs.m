function [ Bwake, Cqwake ] = WakeSourceCoeffs( wake,prf )
%WAKESOURCECOEFFS   calculates 
%                      Bwake=Bij,  i=1,..,N j=N+1,..,N+NW
%                           -> influence of wake sources on airfoil node Eqs
%                      Cqwake=Cij,  i=N+1,..,N+NW j=N+1,..,N+NW
%                           -> influence of wake sources on wake node Eqs


e=prf.panels.e(1:end,:);
nw=wake.n;
ew=wake.e;


NW=evalin('base','NW');%length(wake.x);
N=evalin('base','N');

% all nodes on wake and airfoil
xn=[transpose(prf.nodes.X), transpose(prf.nodes.Y)];
xn=[xn; wake.x, wake.y];

% Panel node 1 and 2 coordinates on wake
xp1=wake.x(1:end-1,1);
yp1=wake.y(1:end-1,1);
xp2=wake.x(2:end,1);
yp2=wake.y(2:end,1);


%local KOS
% line index i -> equation, Loadingpoint in node i=1,..,N+NW
% column index j -> panel index j=N+1,..,N+NW  
X=zeros(N+NW,NW-1); % relativ x coordinates from panel start point to Loading point
X2=zeros(N+NW,NW-1); % relativ x coordinates from panel end point to Loading point
Y=zeros(N+NW,NW-1); % relativ y coordinates from panel to Loading point


for i=1:NW+N % sum over all equations
    xi=[xn(i,1)*ones(NW-1,1), xn(i,2)*ones(NW-1,1)];
    r1= xi-[xp1, yp1]; %relativ vector start point - Loading point
    r2= xi-[xp2, yp2]; %relativ vector end point - Loading point
    tmp=transpose(sum(ew.*r1,2)); 
    X(i,:)=tmp;
    tmp=transpose(sum(ew.*r2,2));
    X2(i,:)=tmp;
    tmp=transpose(sum(nw.*r1,2));
    Y(i,:)=tmp;
end

r1=X.^2+Y.^2; % r^2
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

S1= find(r1<1e-10);
lnr1(S1)=0;
S2= find(r2<1e-10);
lnr2(S2)=0;




% airfoil velocity
%---------------------------
sgn1=sign(X(1:1:N,:));
sgn2=sign(X2(1:1:N,:));
ts1=atan2(sgn1.*X(1:1:N,:),sgn1.*Y(1:1:N,:))+pi*(ones(N,NW-1)-sgn1)/2; 
ts2=atan2(sgn2.*X2(1:1:N,:),sgn2.*Y(1:1:N,:))+pi*(ones(N,NW-1)-sgn2)/2;

% correction for angles 
cor=InvAngle(atan2(-e(:,1),e(:,2)) + pi*ones(N,1) ) ;
cor=cor*ones(1,NW-1);
ts1=ts1-cor ;
ts2=ts2-cor ;

% constant ansatz -> only q_j appears in integral for j-th panel
psi = (1/(2*pi))*( -X(1:1:N,:).*ts1 +X2(1:1:N,:).*ts2 + Y(1:1:N,:).*(lnr1(1:1:N,:)-lnr2(1:1:N,:)) );
Bwake=[psi, zeros(N,1)]; % no source contribution on last wake node



% Wake velocity
%---------------------------
k= ew(:,1).*nw(:,2) - ew(:,2).*nw(:,1);
k=transpose(k*ones(1,NW));

Cqwake= k.*(lnr1(N+1:1:N+NW,:)-lnr2(N+1:1:N+NW,:))/(4*pi);
Cqwake=[Cqwake, zeros(NW,1)];%last value doesn´t appear for linear ansatz

% Cij= 0.5*(e1n2-n1t2)*L^2 - X1*(e1n2-n1t2)*L -2*n1*n2*Y*L


end

