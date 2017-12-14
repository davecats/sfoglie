function [ psi ] = evaluateInviscFieldSol(xi,fld,prf,withoutTE)
%EVALUATEINVISCFIELDSOL     evaluates the streamfunktion for an arbitrary
%                           Field point
%                           if withoutTE is comitted, gamma_TE and q_TE will be neglected 


%---------------------------------------------------
n=prf.panels.n(:,:); % normal vector of each panel
e=prf.panels.e(:,:); % tangent vector of each panel
L=prf.panels.L(:,:); % panel length
N=length(L)-1; % number of panels (without TE panel), N+1 number of nodes

xp1=prf.panels.X(:,1);% 1- start point of panel
yp1=prf.panels.Y(:,1);
xp2=prf.panels.X(:,2);% 2- end point of panel
yp2=prf.panels.Y(:,2);





xiv=[xi(1)*ones(N+1,1), xi(2)*ones(N+1,1)]; % Loading point
r1= xiv-[xp1, yp1]; %relativ vector start point - Loading point
r2= xiv-[xp2, yp2]; %relativ vector end point - Loading point   
X=  transpose(sum(e.*r1,2)); % scalarproduct of relativ vektor with panel tangential vector
X2=transpose(sum(e.*r2,2));
Y= transpose(sum(n.*r1,2));  % scalarproduct of relativ vektor with panel normal vector

r1=X.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

% angles for vorticity coeffs
sgn=sign(Y);
t1=atan2(sgn.*X,sgn.*Y)+pi*(ones(1,N+1)-sgn)/2; 
t2=atan2(sgn.*X2,sgn.*Y)+pi*(ones(1,N+1)-sgn)/2; 

% Filter singularity points out
S1= find(r1<1e-10);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-10);
lnr2(S2)=0;
t2(S2)=0;


% Vortex Coefficients
%--------------------------------------------------------------
LM=transpose(L);
psi1=-(1./(4*pi*LM)).*( (X.*(X-2*LM)-Y.^2).*lnr1 - ((LM-X).^2-Y.^2).*lnr2 ...
                 + 2*Y.*(LM-X).*(t2-t1) +LM.*(3*LM/2-X) );                
psi2=-(1./(4*pi*LM)).*( (Y.^2-X.^2).*lnr1 +(X.^2-Y.^2-LM.^2).*lnr2 ...
                 + 2*Y.*X.*(t2-t1) +(1/2)*LM.*(LM+2*X));

% Trailing edge coefficient for gamma_TE              
TEg=psi1(:,end)+psi2(:,end);

% psi1 is coefficient for j-th node in j-th panel -> last node doesn´t appear 
psi1=[psi1(1:end-1), 0];
% psi2 is coefficient for j+1-th node in j-th panel -> first node doesn´t appear 
psi2=[0, psi2(1:end-1)];

% Coefficient matrix for gammas
A=psi1+psi2;

if nargin==3
%insert TE contribution
M=N+1;
TEq=linesource([prf.panels.X(M,1) prf.panels.Y(M,1);  ...
                             prf.panels.X(M,2) prf.panels.Y(M,2)], ... 
                            [xi(1), xi(2)]);


TE=prf.SdotT*TEg-prf.ScrossT*TEq;
A(1)=A(1) + TE;
end

%---------------------------------------------------

psi=A*fld.gamma;
u=evalin('base','ui');
v= evalin('base','vi');
% add freee stream
psi=psi+u*xi(2)-v*xi(1);

% Old version with loop -> slower

% M = length(fld.gamma);    % number of panels (including trailing edge)
% N = M-1;  
% 
% A=zeros(1,N+1);
% for j=1:N 
%      % Coefficients of gamma_j and gamma_j+1 from j-th panel -> linear ansatz 
%      [psi1,psi2]=linevortex([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
%                              prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
%                             [xi(1), xi(2)]); 
%       A([j j+1]) = A([j j+1]) + [psi1 psi2];       
% end
% 
% 
% if nargin==3
% %insert TE contribution
% [psiT1,psiT2]=linevortex([prf.panels.X(M,1) prf.panels.Y(M,1);  ...
%                              prf.panels.X(M,2) prf.panels.Y(M,2)], ... 
%                             [xi(1), xi(2)]); 
% psiTE= psiT1+psiT2;
% A(1)=A(1)+prf.SdotT*psiTE;
% psiS= linesource([prf.panels.X(M,1) prf.panels.Y(M,1);  ...
%                              prf.panels.X(M,2) prf.panels.Y(M,2)], ... 
%                             [xi(1), xi(2)]);
% A(1)=A(1)-prf.ScrossT*psiS ;
% end
% %---------------------------------------------------
% 
% psi=A*fld.gamma;
% u=evalin('base','ui');
% v= evalin('base','vi');
% % add freee stream
% psi=psi+u*xi(2)-v*xi(1);


end

