function [ Cg, Cq ] = GradPsi( prf,wake )
%GRADPSIG  calculates the coefficients for velocity U=grad(Psi) \cdot n for
%          i=1,..,Nw ; j=1,..,N -> airfoil influence on wake node equations
%           -> U_i= Cg_ij gamma_j + Cq_ij q_j


n=prf.panels.n(1:end,:); 
e=prf.panels.e(1:end,:);
L=prf.panels.L(1:end,:);


NW=evalin('base','NW');%length(wake.x);
N=evalin('base','N');

% nodes for collocation method
xn=[wake.x, wake.y];

% Panel node 1 and 2 coordinates 
xp1=prf.panels.X(1:end,1);
yp1=prf.panels.Y(1:end,1);
xp2=prf.panels.X(1:end,2);
yp2=prf.panels.Y(1:end,2);

%local KOS
% line index i -> equation, Loadingpoint in node i=1,..,N+1
% column index j -> panel index j=1,..,N+1  N+1 TE panel
X=zeros(NW,N); % relativ x coordinates from panel start point to Loading point
X2=zeros(NW,N); % relativ x coordinates from panel end point to Loading point
Y=zeros(NW,N); % relativ y coordinates from panel to Loading point

for i=1:NW % sum over wake equations
    xi=[xn(i,1)*ones(N,1), xn(i,2)*ones(N,1)];
    r1= xi-[xp1, yp1]; %relativ vector start point - Loading point
    r2= xi-[xp2, yp2]; %relativ vector end point - Loading point
    tmp=transpose(sum(e.*r1,2)); 
    X(i,:)=tmp;
    tmp=transpose(sum(e.*r2,2));
    X2(i,:)=tmp;
    tmp=transpose(sum(n.*r1,2));
    Y(i,:)=tmp;
end

r1=X.^2+Y.^2; % r^2
r2=X2.^2+Y.^2;
%Loadingpoint only for wake -> no singularity here
lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;
sgn=sign(Y);
t1=atan2(sgn.*X,sgn.*Y)+pi*(ones(NW,N)-sgn)/2; 
t2=atan2(sgn.*X2,sgn.*Y)+pi*(ones(NW,N)-sgn)/2; 

absn=sum(n.*n,2); % |n|^2
absn= transpose(absn*ones(1,NW));

% circulation
%---------------------------
LM=transpose(L*ones(1,NW));
A=1./(4*pi*LM).*( absn.*Y.*(lnr1-lnr2)+2*absn.*(X-LM).*(t2-t1) );
GTE=A(:,end);% extract trailing edge panel
A=[A(:,1:end-1), ones(NW,1)]; % coefficient for j-th node in j-th panel -> last node doesn´t appear 
B=1./(4*pi*LM).*( Y.*(lnr1-lnr2)+2*X.*absn.*(t2-t1) );
GTE=(GTE+B(:,end));% extract trailing edge panel
B=[ones(NW,1),B(:,1:end-1)]; % coefficient for j+1-th node in j-th panel -> first node doesn´t appear 

Cg=A+B;
% a=(Y|n|^2(lnr1-lnr2) +2(X-L)|n|^2 (t2-t1) )/(4 pi L)
% b=(Y(lnr1-lnr2) +2X|n|^2 (t2-t1) )/(4 pi L)
% Integral from j-th panel:  a_j gamma_j +b_j gamma_j+1
% Indexschift -> Cij= (a_j + b_j-1)   j=2,...,N-1
%             -> Ci1= a_1  
%             -> CiN= b_N-1  


% source
%---------------------------
k= e(:,1).*n(:,2) - e(:,2).*n(:,1);
k=transpose(k*ones(1,NW));


Cq= k.*(lnr1-lnr2)/(4*pi);
% Cij= (e1n2-n1t2)(lnr1-lnr2)/(4*pi)

STE=Cq(:,end);
Cq= [Cq(:,1:end-1), zeros(NW,1) ]; %last value doesn´t appear for linear ansatz

% Add Trailing edge influence to matrix A
TE= prf.SdotT*GTE-prf.ScrossT*STE;
Cg(:,1)=Cg(:,1)+TE;
    
end

