
close all
clear all

addpath('./panel/')
addpath('./Amplification/')
addpath('./geometry/')
addpath('./utilities/')
addpath('./models/')
addpath('./wake/')
addpath('./BoundaryLayer/')
set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');

% sets beta values to get solutions
beta= [-0.1988, -0.19:0.01:2 ];
%beta= [-0.1988, -0.198:0.001:2 ];
H12=zeros(size(beta));
H32=H12;Cf=H12;CD2=H12; ddf0=H12;
b1=H12;b2=H12;b3=H12;bD=H12;

for i=1:length(beta)
    [~,erg ]=FalknerSkanInt(beta(i),false,false,false);
    H12(i)=erg.H12;
    H32(i)=erg.H32;
    Cf(i) =erg.Cf;
    CD2(i)=2*erg.CD;
    b1(i) = erg.beta1;
    b2(i) = erg.beta2;
    b3(i) = erg.beta3;
    bD(i) = erg.betaD;
    ddf0(i)=erg.s;
end


CD=CD2/2;
m= beta./(2-beta);
% FSerg= [beta;m; H12; H32; Cf; CD; ddf0; b1;b2;b3;bD];
% dlmwrite('Falkner-Skan-solutions.txt',FSerg);

 
% Plots
%%

% correlations from Falkner-Skan-solution
figure
hold on
plot(beta,H12)
plot(beta,H32)
plot(beta,Cf)
plot(beta,CD2)
legend('$H_{12}$','$H_{32}$ ','$C_f$','$C_D$','location','best');   


figure
hold on
plot(H12,H32)
plot(H12,beta)
plot(H12,Cf)  
plot(H12,CD2)
legend('$H_{32}$','$\beta$ ','$C_f$','$C_D$','location','best');   

% comparisson with models from Drela 1989

CFT=2*CF2lam( H12,ones(size(H12)) );
HST=H32lam(H12);
CDT=CD2lam( H12,ones(size(H12)) ).*HST;
mT= (0.058*(H12-4.0).^2./(H12-1.0) - 0.068) ./ (6.54*H12 - 14.07) .* H12.^2;
BT= 2*mT./(mT+1);
% figure
% hold on
% plot(H12,HST)
% plot(H12,BT)
% plot(H12,CFT)
% plot(H12,CDT)
% legend('$H_{32}$','$\beta$ ','$C_f$','$C_D$','location','best');  



% comparisson with older models from Drela 1987


CFT2=2*CF2lam( H12,ones(size(H12)) ,'old');
HST2=H32lam(H12,'old');
CDT2=CD2lam( H12,ones(size(H12)),'old' ).*HST2;
% figure
% hold on
% plot(H12,HST2)
% plot(H12,BT)
% plot(H12,CFT2)
% plot(H12,CDT2)
% legend('$H_{32}$','$\beta$ ','$C_f$','$C_D$','location','best');  


% comparisson with models from Eppler 1980

ind=find(H32>=1.51509 & H32<1.57258);
ind2=find(H32>=1.57258);

HE=79.870845 - 89.582142*H32 + 25.715786*H32.^2;
HE(ind)= 4.02922 - (583.60182 - 724.55916*H32(ind) + 227.1822*H32(ind).^2).*sqrt(H32(ind)-1.51509);

CfE=1.372391 - 4.226253*H32 + 2.221687*H32.^2;
CfE(ind)= 2.512589 - 1.686095*HE(ind) + 0.391541*HE(ind).^2 - 0.031729*HE(ind).^3;

CfE=2*CfE;

CDE= 2* (7.853976 - 10.260551*H32 + 3.418898* H32.^2 );

% figure
% hold on
% plot(HE,H32)
% plot(HE,beta)
% plot(HE,CfE)  
% plot(HE,CDE)
% legend('$H_{32}$','$\beta$ ','$C_f$','$C_D$','location','best');   

% compare
%-------------------------------------------------------
figure
hold on
title('Vergleich: H_3_2');
plot(H12,H32)
plot(H12,HST)
plot(H12,HST2)
plot(HE,H32)
legend('Falkner skan sol','Drela 1989 ','Drela 1987','Eppler 1980','location','best');  


figure
hold on
title('Vergleich: C_f');
plot(H12,Cf)
plot(H12,CFT)
plot(H12,CFT2)
plot(HE,CfE)
legend('Falkner skan sol','Drela 1989 ','Drela 1987','Eppler 1980','location','best');  

figure
hold on
plot(H12,CD2)
plot(H12,CDT)
plot(H12,CDT2)
plot(HE,CDE)
title('Vergleich: 2 C_D');
legend('Falkner skan sol','Drela 1989 ','Drela 1987','Eppler 1980','location','best');  











