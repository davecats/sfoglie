function [dpsi_dxi,dpsi_deta] = GradPsi(xi,eta, prf )
%GRADPSI 




L=prf.panels.L; 
N=length(xi);


% node coordinates 
[X1,X2,Y]=GetLocalRelCoord(xi,eta,prf);


r1=X1.^2 + Y.^2; % r^2 
r2=X2.^2 + Y.^2;

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


LM=ones(N,1)*L;

% partial derivates
dps1_dX1=   2*(X1-LM).*(lnr1-lnr2) + 2*Y.*(t1-t2)-LM  ;
dps1_dY =   2*Y.*(lnr2-lnr1) + 2*(LM-X1).*(t2-t1);
             

dps1_dlnr1= (X1.*(X1-2*LM) - Y.^2);
dps1_dlnr2= - ( (LM-X1).^2 - Y.^2 );

dps1_dt1= -2*Y.*(LM-X1);
dps1_dt2=  2*Y.*(LM-X1);


dps2_dX1=   2*X1.*(lnr2-lnr1) + 2*Y.*(t2-t1) + LM  ;
dps2_dY =   2*Y.*(lnr1-lnr2) + 2*X1.*(t2-t1) ;


dps2_dlnr1= (X1.*(X1-2*LM) - Y.^2);
dps2_dlnr2= - ( (LM-X1).^2 - Y.^2 );

dps2_dt1= -2*Y.*(LM-X1);
dps2_dt2=  2*Y.*(LM-X1);


dlnr1_dX = 2*X1./r1;
dlnr2_dX = 2*X2./r2;
dlnr1_dY = 2*Y ./r1;
dlnr2_dY = 2*Y ./r2;

dt1_dX =   Y./r1;
dt1_dY = -X1./r1;

dt2_dX =   Y./r2;
dt2_dY = -X2./r2;           
           
% total derivates
dps1_dX1= -(1./(4*pi*LM)).*(dps1_dX1 + dps1_dlnr1.*dlnr1_dX + dps1_dlnr2.*dlnr2_dX + dps1_dt1.*dt1_dX + dps1_dt2.*dt2_dX);
dps1_dY = -(1./(4*pi*LM)).*(dps1_dY  + dps1_dlnr1.*dlnr1_dY + dps1_dlnr2.*dlnr2_dY + dps1_dt1.*dt1_dY + dps1_dt2.*dt2_dY);

dps2_dX1= -(1./(4*pi*LM)).*(dps2_dX1 + dps2_dlnr1.*dlnr1_dX + dps2_dlnr2.*dlnr2_dX + dps2_dt1.*dt1_dX + dps2_dt2.*dt2_dX);
dps2_dY = -(1./(4*pi*LM)).*(dps2_dY  + dps2_dlnr1.*dlnr1_dY + dps2_dlnr2.*dlnr2_dY + dps2_dt1.*dt1_dY + dps2_dt2.*dt2_dY);

dX1_dxi = ones(N,1)*prf.panels.e(1,:);
dX1_deta= ones(N,1)*prf.panels.e(2,:);

dY_dxi = ones(N,1)*prf.panels.n(1,:);
dY_deta= ones(N,1)*prf.panels.n(2,:);

dps1_dxi = dps1_dX1.*dX1_dxi  + dps1_dY.*dY_dxi;
dps1_deta= dps1_dX1.*dX1_deta + dps1_dY.*dY_deta;

dps2_dxi = dps2_dX1.*dX1_dxi  + dps2_dY.*dY_dxi;
dps2_deta= dps2_dX1.*dX1_deta + dps2_dY.*dY_deta;

dGTE_dxi =(dps1_dxi(:,end) + dps2_dxi(:,end) );
dGTE_deta=(dps1_deta(:,end)+ dps2_deta(:,end));


dQTE_dX1= (1/(2*pi))*( t2(:,end) - t1(:,end) - X1(:,end).*dt1_dX(:,end) + X2(:,end).*dt2_dX(:,end) +Y(:,end).*(dlnr1_dX(:,end) - dlnr2_dX(:,end)) );
dQTE_dY = (1/(2*pi))*( X2(:,end).*dt2_dY(:,end) - X1(:,end).*dt1_dY(:,end) +Y(:,end).*(dlnr1_dY(:,end) - dlnr2_dY(:,end)) + (lnr1(:,end) - lnr2(:,end)) );

dQTE_dxi = dQTE_dX1.*dX1_dxi(:,end)  + dQTE_dY.*dY_dxi(:,end); 
dQTE_deta= dQTE_dX1.*dX1_deta(:,end) + dQTE_dY.*dY_deta(:,end); 

p1_xi=[dps1_dxi(:,1:end-1), zeros(N,1)];
p2_xi=[zeros(N,1), dps2_dxi(:,1:end-1)];
p1_eta=[dps1_deta(:,1:end-1), zeros(N,1)];
p2_eta=[zeros(N,1), dps2_deta(:,1:end-1)];

dpsi_dxi = p1_xi + p2_xi;
dpsi_deta= p1_eta + p2_eta;

if prf.sharpTE
    dGTE_dxi=0;
    dGTE_deta=0;
end

TE_xi=  prf.SdotT*dGTE_dxi + prf.ScrossT*dQTE_dxi;
TE_eta=  prf.SdotT*dGTE_deta + prf.ScrossT*dQTE_deta;

dpsi_dxi(:,1)  =dpsi_dxi(:,1)   + TE_xi/2;
dpsi_dxi(:,end)=dpsi_dxi(:,end) - TE_xi/2;

dpsi_deta(:,1)  =dpsi_deta(:,1)   + TE_eta/2;
dpsi_deta(:,end)=dpsi_deta(:,end) - TE_eta/2;



end

