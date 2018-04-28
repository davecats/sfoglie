function [ psi ] = evaluateInviscFieldSol(xi,fld,prf,withoutTE)
%EVALUATEINVISCFIELDSOL     evaluates the streamfunktion for an arbitrary
%                           Field point
%                           if withoutTE is comitted, gamma_TE and q_TE will be neglected 


%---------------------------------------------------

if nargin==3
    withoutTE=false;
end
    

L=prf.panels.L; % panel length

% node coordinates 
[X,X2,Y]=GetLocalRelCoord(xi(1),xi(2),prf);

r1=X.^2  + Y.^2; % r^2 
r2=X2.^2 + Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;


Y=-Y;
t1=atan2(X ,Y); 
t2=atan2(X2,Y); 

% Filter singularity points out
S1= find(r1<1e-10);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-10);
lnr2(S2)=0;
t2(S2)=0;


% Vortex Coefficients
%--------------------------------------------------------------

psi1=-(1./(4*pi*L)).*( (X.*(X-2*L)-Y.^2).*lnr1 - ((L-X).^2-Y.^2).*lnr2 ...
                 + 2*Y.*(L-X).*(t2-t1) +L.*(3*L/2-X) );                
psi2=-(1./(4*pi*L)).*( (Y.^2-X.^2).*lnr1 +(X.^2-Y.^2-L.^2).*lnr2 ...
                 + 2*Y.*X.*(t2-t1) +(1/2)*L.*(L+2*X));

if prf.sharpTE==false && withoutTE==false;
    % Trailing edge coefficient for gamma_TE              
    GTE=psi1(end)+psi2(end);

elseif prf.sharpTE
    GTE=(1/(2*pi))*Y(:,end)*pi;
    
end
% psi1 is coefficient for j-th node in j-th panel -> last node doesn´t appear 
psi1=[psi1(1:end-1), 0];
% psi2 is coefficient for j+1-th node in j-th panel -> first node doesn´t appear 
psi2=[0, psi2(1:end-1)];
    
    
% Coefficient matrix for gammas
A=psi1+psi2;

if nargin==3 && prf.sharpTE==false
    %insert TE contribution
    cor=prf.panels.theta(end);
    t1=t1(end)-cor;
    t2=t2(end)-cor;

    qTE = (1/(2*pi))*( - X(end).*t1 + X2(end).*t2 + Y(end).*(lnr1(end) - lnr2(end)) );

    TE=prf.SdotT*GTE + prf.ScrossT*qTE;
    A(1)=A(1) + TE;
    
end

%---------------------------------------------------

psi=A*fld.gamma; 


% add freee stream
psi=psi + fld.ui*xi(2)-fld.vi*xi(1);


end

