function [ HS,dHS_dH,dHS_dRet  ] = H32turb( H, Ret )
%H32TURB calculates H32 and its the derivate in respect to the shape parameter H12 for turbulent case


% Help parameter H0
dim=size(H);
n=length(H);
H0=4*ones(dim);
dH0_dRet=zeros(dim);
H0(Ret>400)=3 + 400/Ret(Ret>400);
dH0_dRet(Ret>400)=-400/Ret(Ret>400).^2;


% threatment of small Rets
Ret(Ret<200)=200;
dRet_dRet=ones(n,1);
dRet_dRet(Ret<200)=0;


HS= zeros(dim);
dHS_dH  = zeros(dim);
dHS_dRet= zeros(dim);

% H12<H0 -> attached branch
ind1=find(H<=H0);

HR=(H0(ind1)-H(ind1))./(H0(ind1)-1);
dHR_dH=-1./(H0(ind1)-1);
dHR_dRet=(1-HR)./(H0(ind1)-1).*dH0_dRet(ind1);

HS(ind1)=1.5 + 4./Ret(ind1) + ( 0.5-4./Ret(ind1) ).*HR.^2 * 1.5./(H(ind1)+0.5);

dHS_dH(ind1)=   - ( 0.5-4./Ret(ind1) ).*HR.^2 .* 1.5./(H(ind1)+0.5).^2 ...
                + ( 1  -8./Ret(ind1) ).*HR    .* 1.5./(H(ind1)+0.5) .*dHR_dH;
            
dHS_dRet(ind1)=  (HR.^2*1.5./(H(ind1)+0.5) -1).*4./Ret(ind1).^2.*dRet_dRet(ind1) ...
                + ( 1  -8./Ret(ind1) ).*HR    .* 1.5./(H(ind1)+0.5) .*dHR_dRet;            
            
% H12>H0 -> seperation
ind2=find(H>H0);


logRet=log(Ret(ind2));
delta = H(ind2)-H0(ind2);
k1 = delta + 4./logRet;
k2 = 0.007*logRet./k1.^2 + 0.015./H(ind2);
dk2_dH  =-0.014*logRet./k1.^3 - 0.015./H(ind2).^2;
dk2_dRet=-0.014*logRet./k1.^3 .* (-dH0_dRet(ind2) - 4./(logRet.^2.*Ret(ind2)).*dRet_dRet(ind2) ) ...
            + 0.007./(k1.^2.*Ret(ind2)).*dRet_dRet(ind2) ;

HS(ind2)= delta.^2 .* k2 + 1.5 + 4./Ret(ind2);

dHS_dH(ind2)= 2*delta.*k2 + delta.^2.*dk2_dH;

dHS_dRet(ind2)= delta.^2.*dk2_dRet - 4./Ret(ind2).^2.*dRet_dRet(ind2) ...
                - 2*delta.*k2.*dH0_dRet(ind2);

                
end

