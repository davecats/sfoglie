function [ A,B] = VortexAndSourceDistr( prf )
%VORTEXANDSOURCEDISTR  calculates the Matrix of source and vorticity
%                      contributions for the collocation method
%                      line index i: Index of equation -> Loading point at i-th node
%                      column index j: Index for influence of j-th node
%                                      value to boundary integral


% with TE panel
n=prf.panels.n(1:end,:); % normal vector of each panel
e=prf.panels.e(1:end,:); % tangent vector of each panel
L=prf.panels.L(1:end,:); % panel length
N=length(L)-1; % number of panels (without TE panel), N+1 number of nodes

% node coordinates
xn=[transpose(prf.nodes.X), transpose(prf.nodes.Y)];

% Panel node 1 and 2 coordinates
xp1=prf.panels.X(:,1);% 1- start point of panel
yp1=prf.panels.Y(:,1);
xp2=prf.panels.X(:,2);% 2- end point of panel
yp2=prf.panels.Y(:,2);

%local KOS
% line index i -> equation, Loadingpoint in node i=1,..,N+1
% column index j -> panel index j=1,..,N+1 - N+1 TE panel
X=zeros(N+1,N+1); % relativ x coordinates from panel start point to Loading point
X2=zeros(N+1,N+1); % relativ x coordinates from panel end point to Loading point
Y=zeros(N+1,N+1); % relativ y coordinates from panel to Loading point


for i=1:N+1 % sum over all equations
    xi=[xn(i,1)*ones(N+1,1), xn(i,2)*ones(N+1,1)]; % Loading point
    r1= xi-[xp1, yp1]; %relativ vector start point - Loading point
    r2= xi-[xp2, yp2]; %relativ vector end point - Loading point
    tmp=transpose(sum(e.*r1,2));    % scalarproduct of relativ vektor with 
    X(i,:)=tmp;                     % panel tangential vector
    tmp=transpose(sum(e.*r2,2));
    X2(i,:)=tmp;
    tmp=transpose(sum(n.*r1,2));    % scalarproduct of relativ vektor with 
    Y(i,:)=tmp;                     % panel normal vector
end

r1=X.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

% angles for vorticity coeffs
sgn=sign(Y);
t1=atan2(sgn.*X,sgn.*Y)+pi*(ones(N+1,N+1)-sgn)/2; 
t2=atan2(sgn.*X2,sgn.*Y)+pi*(ones(N+1,N+1)-sgn)/2; 

% Filter singularity points out
S1= find(r1<1e-10);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-10);
lnr2(S2)=0;
t2(S2)=0;


% Vortex Coefficients
%--------------------------------------------------------------
LM=transpose(L*ones(1,N+1));
psi1=-(1./(4*pi*LM)).*( (X.*(X-2*LM)-Y.^2).*lnr1 - ((LM-X).^2-Y.^2).*lnr2 ...
                 + 2*Y.*(LM-X).*(t2-t1) +LM.*(3*LM/2-X) );                
psi2=-(1./(4*pi*LM)).*( (Y.^2-X.^2).*lnr1 +(X.^2-Y.^2-LM.^2).*lnr2 ...
                 + 2*Y.*X.*(t2-t1) +(1/2)*LM.*(LM+2*X));

% Trailing edge coefficient for gamma_TE              
GTE=psi1(:,end)+psi2(:,end);

% psi1 is coefficient for j-th node in j-th panel -> last node doesn´t appear 
psi1=[psi1(:,1:end-1), zeros(N+1,1)];
% psi2 is coefficient for j+1-th node in j-th panel -> first node doesn´t appear 
psi2=[zeros(N+1,1), psi2(:,1:end-1)];

% Coefficient matrix for gammas
A=psi1+psi2;

% source Coefficients
%--------------------------------------------------------------

% angle for source coeffs
sgn1=sign(X);
sgn2=sign(X2);
ts1=atan2(sgn1.*X,sgn1.*Y)+pi*(ones(N+1,N+1)-sgn1)/2; 
ts2=atan2(sgn2.*X2,sgn2.*Y)+pi*(ones(N+1,N+1)-sgn2)/2;
ts1(S1)=0;ts2(S2)=0;
% correction for angles 
cor=InvAngle(atan2(-e(:,1),e(:,2)) + pi*ones(N+1,1) ) ;
cor=transpose(cor*ones(1,N+1));
ts1=ts1-cor ;
ts2=ts2-cor ;

% constant ansatz -> only q_j appears in integral for j-th panel
psi = (1/(2*pi))*( -X.*ts1 +X2.*ts2 + Y.*(lnr1-lnr2) );
B=psi;


% Trailing edge coefficient for q_TE  
STE=psi(:,end);


% Add Trailing edge influence to matrix A

TE= prf.SdotT*GTE-prf.ScrossT*STE; % s_t gamma_TE + scrosst q_TE
A(:,1)=A(:,1)+TE; % kutta Condition already applied here


end

