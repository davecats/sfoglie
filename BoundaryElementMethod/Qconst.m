function B = Qconst( prf, wake)
%QCONST  calculates source coefficients B with constant ansatz
%       if only prf is comitted -> Bij, i=1,..,N ; ij=1,..,N
%       if wake is comitted     -> Bij, i=1,..,N ; ij=N+1,..,N+NW


if nargin==1 % Contribution of airfoil nodes
    L=prf.panels.L(1:end-1); % panel length without TE panel
    N=length(L); % number of panels (without TE panel), N+1 number of nodes


    % panel angle
    panAng=prf.panels.theta(1:end-1);

    % node coordinates 
    [X1,X2,Y]=GetLocalRelCoord(transpose(prf.nodes.X),transpose(prf.nodes.Y),prf,'noTE');
else % Contribution of wake nodes
    L=wake.L;
    N=length(prf.nodes.X)-1;
    panAng=wake.theta;

    [X1,X2,Y]=WakeRelCoord(transpose(prf.nodes.X),transpose(prf.nodes.Y),wake);
end



r1=X1.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;


lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;


t1=atan2(X1,Y); 
t2=atan2(X2,Y); 


% t1=atan(X1./Y); 
% t2=atan(X2./Y); 

% sgn=ones(N+1,N+1);
% sgn(1,2:end)=sign(Y(1,2:end));
% sgn(N+1,1:end-1)=sign(Y(N+1,1:end-1));
% t1=atan2(sgn.*X1,sgn.*Y) + (ones(N+1,N+1) - sgn)*pi/2;
% t2=atan2(sgn.*X2,sgn.*Y) + (ones(N+1,N+1) - sgn)*pi/2;


% Filter singularity points out
S1= find(r1<1e-11);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-11);
lnr2(S2)=0;
t2(S2)=0;

cor=ones(N+1,1)*panAng;
t1=t1-cor; 
t2=t2-cor; 

B= (1/(2*pi))*( -X1.*t1 + X2.*t2 + Y.*(lnr1 - lnr2) );
B= [B, zeros(size(B(:,1)))];



end

