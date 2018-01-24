function [ A,B] = VortexAndSourceDistr( prf )
%VORTEXANDSOURCEDISTR  calculates the Matrix of source and vorticity
%                      contributions for the collocation method
%                      line index i: Index of equation -> Loading point at i-th node
%                      column index j: Index for influence of j-th node
%                                      value to boundary integral


% with TE panel
L=prf.panels.L; % panel length
if prf.IsSharp
    N=length(L);
else
    N=length(L)-1; % number of panels (without TE panel), N+1 number of nodes
end

 % node coordinates 
 [X,X2,Y]=GetLocalRelCoord(transpose(prf.nodes.X),transpose(prf.nodes.Y),prf);


r1=X.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

t1=atan2(X ,Y);
t2=atan2(X2,Y);

% Filter singularity points out
S1= find(r1<1e-11);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-11);
lnr2(S2)=0;
t2(S2)=0;

% t1X=load('./XFoilWerte/t1.txt');
% t2X=load('./XFoilWerte/t2.txt');


% Vortex Coefficients
%--------------------------------------------------------------
LM=ones(N+1,1)*L;
psi1=-(1./(4*pi*LM)).*( (X.*(X-2*LM)-Y.^2).*lnr1 - ((LM-X).^2-Y.^2).*lnr2 ...
                 + 2*Y.*(LM-X).*(t2-t1) +LM.*(3*LM/2-X) );                
psi2=-(1./(4*pi*LM)).*( (Y.^2-X.^2).*lnr1 +(X.^2-Y.^2-LM.^2).*lnr2 ...
                 + 2*Y.*X.*(t2-t1) +(1/2)*LM.*(LM+2*X));

             
% Trailing edge coefficient for gamma_TE                  
GTE=(psi1(:,end)+ psi2(:,end));
% psi1 is coefficient for j-th node in j-th panel -> last node doesn´t appear 
psi1=[psi1(:,1:end-1), zeros(N+1,1)];
% psi2 is coefficient for j+1-th node in j-th panel -> first node doesn´t appear 
psi2=[zeros(N+1,1), psi2(:,1:end-1)];

% Coefficient matrix for gammas
A=psi1+psi2;

% source Coefficients
%--------------------------------------------------------------

% correction with panel angles
cor=ones(N+1,1)*prf.panels.theta;
t1=t1-cor;
t2=t2-cor;

% without TE panel
% constant ansatz -> only q_j appears in integral for j-th panel
psi = (1/(2*pi))*( -X.*t1 + X2.*t2 + Y.*(lnr1 - lnr2) );
B=psi(:,1:end-1);
qTE=psi(:,end);


    
%     % TE circulation
%      GTE =  (1/(2*pi))*( X(:,end).*lnr1(:,end) -  X2(:,end).*lnr2(:,end) ...
%                          + Y(:,end).*(t2(:,end) - t1(:,end)) - LM(:,end) );
%     % TE source
%     qTEt = (1/(2*pi))*( - X(:,end).*t1(:,end) + X2(:,end).*t2(:,end) ...
%                        + Y(:,end).*(lnr1(:,end) - lnr2(:,end)) );

    
   
    % Add Trailing edge influence to matrix A
    TE=  prf.SdotT*GTE + prf.ScrossT*qTE; % s_t gamma_TE + scrosst q_TE
    %A(:,1)=A(:,1)+TE; % kutta Condition already applied here

    A(:,1)=A(:,1)+TE/2; 
    A(:,end)=A(:,end)-TE/2;





    
% % angle for source coeffs
% Y=-Y;
% sgn1=sign(X);
% sgn2=sign(X2);
% ts1=atan2(sgn1.*X,sgn1.*Y)+pi*(ones(N+1,1)-sgn1)/2; 
% ts2=atan2(sgn2.*X2,sgn2.*Y)+pi*(ones(N+1,1)-sgn2)/2;
% ts1(r1(:,end)<1e-11)=0;ts2(r2(:,end)<1e-11)=0;
% % correction for angles 
% cor=InvAngle(atan2(-e(1,:),e(2,:)) + pi*ones(1,N+1) ) ;
% cor=ones(N+1,1)*cor;
% ts1=ts1-cor(:,end) ;
% ts2=ts2-cor(:,end) ;
% % constant ansatz -> only q_j appears in integral for j-th panel
% qTE = (1/(2*pi))*( -X.*ts1 +X2.*ts2 + Y.*(lnr1(:,end)-lnr2(:,end)) ); 
    
    




end

