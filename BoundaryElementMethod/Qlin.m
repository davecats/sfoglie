function M = Qlin( xi, eta, handle,IsWake)
%QLIN  calculates the airfoil panel Coefficients Bij for Loading point xi
%       input:
%       -  xi = xi_i , Loadingpoint x-coordinate
%       -  eta=eta_i , Loadingpoint y-coordinate
%       -> IsWake=false -> handle=prf   i=1,..,N    -> Bij, j=1,..,N (default)
%       -> IsWake=true  -> handle=wake  i=N,..,N    -> Bij, j=1,..,N
%       

if nargin==3
    IsWake=false;
end
    
if ~IsWake
    [X1,X2,Y]=GetLocalRelCoord( xi,eta,handle,true);
    L=handle.panels.L(1:end-1); % panel length without TE panel

    
    % panel angle
    panAng=handle.panels.theta(1:end-1);
else
    [X1,X2,Y]=WakeRelCoord( xi,eta,handle);
    L=handle.L; % panel length without TE panel
    % panel angle
    panAng=handle.theta; 
end

n=length(X1(:,1));
m=length(X1(1,:));

% panel midpoint
XM=(X1+X2)/2;


r1=X1.^2 + Y.^2; % r^2 
r2=X2.^2 + Y.^2;
rM=XM.^2 + Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;
lnrM=log(rM)/2; 

t1=atan2(X1,Y); 
t2=atan2(X2,Y); 
tM=atan2(XM,Y);


% Filter singularity points out
S1= find(r1<1e-11);
lnr1(S1)=0;
t1(S1)=0;
S2= find(r2<1e-11);
lnr2(S2)=0;
t2(S2)=0;

% adds jth panel angle to each column 
cor=ones(n,1)*panAng; 
t1=t1 - cor; 
t2=t2 - cor; 
tM=tM - cor;


% length of j-1 panel + j panel
Lm= [L(1), L(1:end-1) + L(2:end)]; 
hm=ones(n,1)*Lm;
% length of j panel
h=ones(n,1)*L;
% length of j+1 panel + j panel
Lp= [L(1:end-1) + L(2:end), L(end)]; 
hp=ones(n,1)*Lp;


 % Approximation of first distribution
pp= XM.*tM - X1.*t1 + Y.*(lnr1-lnrM);
pm= ( (X1+XM).*pp + r1.*t1 - rM.*tM + Y.*(XM-X1) )./(X1-XM);



Cjm=  ( -pp + pm )./(4*pi*hm);
Cj = -(  pp + pm )./(4*pi*h);
Cjp=  (  pp.*(1./h + 1./hm) + pm.*(1./h - 1./hm) )/(4*pi);


Cj =[Cj ,zeros(n,1)]; 
Cjm=[Cjm(:,1) + Cjm(:,2), Cjm(:,3:end) ,zeros(n,2)];
Cjp=[zeros(n,1),Cjp];

M1= Cjm + Cj + Cjp;

% % Approximation of second distribution
pp= X2.*t2 - XM.*tM + Y.*(lnrM-lnr2);
pm= ( (X2+XM).*pp + rM.*tM - r2.*t2 + Y.*(X2-XM) )./(XM-X2);


Cj  = -(  pp.*(1./h + 1./hp) + pm.*(1./hp - 1./h) )/(4*pi);
Cjp =  (  pp - pm )./(4*pi*h);
Cjpp=  (  pp + pm )./(4*pi*hp);

Cj =[Cj ,zeros(n,1)];
Cjp=[zeros(n,1),Cjp];
Cjpp=[zeros(n,2), Cjpp(:,1:end-2) ,Cjpp(:,end-1)+Cjpp(:,end)];

M2= Cj + Cjp + Cjpp;


M=(M1+M2);

end

