function [ Cqw ] = GradPsiW( wake, ansatz )
%GRADPSIW calculates the coefficients for velocity U=grad(Psi) \cdot n for
%          i=1,..,Nw ; j=N+1,..,N+NW -> wake influence on wake node equations


NW=length(wake.x);

[X1,X2,Y]=WakeRelCoord(wake.x,wake.y,wake);

r1=X1.^2+Y.^2; % r^2 
r2=X2.^2+Y.^2;

lnr1=log(r1)/2; % ln(r) =ln(r^2)/2
lnr2=log(r2)/2;

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
% k2=e1*n2-e2*n1 =  (n1^2 + n2^2)
k2= wake.e(1,:).*wake.n(2,:) -  wake.e(2,:).*wake.n(1,:);
k2=ones(NW,1)*k2;



if nargin==2
    % constant ansatz
    Cqw=   k1.*(t2-t1) + k2.*(lnr1-lnr2)/(2*pi);
    Cqw=[Cqw, zeros(NW,1)];%last value doesn´t appear for constant ansatz
else % piecewise linear ansatz
    XM=(X1+X2)/2;
    rM=XM.^2+Y.^2;
    lnrM=log(rM)/2; 
    tM=atan2(XM,Y);
    N=length(wake.L);
    
    cor=ones(N+1,1)*wake.theta;
    t1=t1-cor;t2=t2-cor;tM=tM-cor;
    
    % length of j-1 panel + j panel
    Lm= [wake.L(1), wake.L(1:end-1) + wake.L(2:end)]; 
    hm=ones(N+1,1)*Lm;
    % length of j panel
    h=ones(N+1,1)*wake.L;
    % length of j+1 panel + j panel
    Lp= [wake.L(1:end-1)+wake.L(2:end), wake.L(end)]; 
    hp=ones(N+1,1)*Lp;
    
    
     % Approximation of first distribution
    pp= XM.*tM - X1.*t1 + Y.*(lnr1-lnrM);
    pm= ( (X1+XM).*pp +r1.*t1 - rM.*tM + Y.*(XM-X1) )./(X1-XM);

    dpm_dX1= (-(X1+XM).*t1 + pp + 2*X1.*t1 - pm)./(X1-XM);
    dpm_dXM= ( (X1+XM).*tM + pp - 2*XM.*tM + pm)./(X1-XM);
    dpm_dY = ( (X1+XM).*(lnr1-lnrM) + 2*(XM-X1 + Y.*(t1-tM)) )./(X1-XM);

    gpp = k1.*(tM-t1) + k2.*(lnr1-lnrM);
    gpm = dpm_dX1.*k1 + dpm_dXM.*k1 + dpm_dY.*k2;
    
    Cjm=  ( -gpp + gpm )./(4*pi*hm);
    Cj = -(  gpp + gpm )./(4*pi*h);
    Cjp=  (  gpp.*(1./h+1./hm) + gpm.*(1./h-1./hm) )/(4*pi);

    Cj =[Cj ,zeros(N+1,1)];
    Cjm=[Cjm(:,1)+Cjm(:,2),Cjm(:,3:end) ,zeros(N+1,2)];
    Cjp=[zeros(N+1,1),Cjp];

    C1= Cjm + Cj + Cjp;
    
    % % Approximation of second distribution
    pp= X2.*t2 - XM.*tM + Y.*(lnrM-lnr2);
    pm= ( (X2+XM).*pp +rM.*tM - r2.*t2 + Y.*(X2-XM) )./(XM-X2);

    dpm_dXM= (-(X2+XM).*tM + pp + 2*XM.*tM - pm)./(XM-X2);
    dpm_dX2= ( (X2+XM).*t2 + pp - 2*X2.*t2 + pm)./(XM-X2);
    dpm_dY =  ((X2+XM).*(lnrM-lnr2) + 2*(X2-XM + Y.*(tM-t2)) )./(XM-X2);

    gpp = k1.*(t2-tM) + k2.*(lnrM-lnr2);
    gpm = dpm_dXM.*k1 + dpm_dX2.*k1 + dpm_dY.*k2;
    
    Cj  = -(  gpp.*(1./h+1./hp) + gpm.*(1./hp-1./h) )/(4*pi);
    Cjp =  (  gpp - gpm )./(4*pi*h);
    Cjpp=  (  gpp + gpm )./(4*pi*hp);

    Cj =[Cj ,zeros(N+1,1)];
    Cjp=[zeros(N+1,1),Cjp];
    Cjpp=[zeros(N+1,2), Cjpp(:,1:end-2) ,Cjpp(:,end-1)+Cjpp(:,end)];

    C2= Cj + Cjp + Cjpp;


    Cqw=(C1+C2);
end







end

