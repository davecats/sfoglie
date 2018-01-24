function [ CD2, dCD2_dH,dCD2_dRet ] = CD2Lwake( H12,Ret )
%CD2LWAKE calculates CD for laminar wake


k=(1.1.*(1-1./H12).^2 )./H12;
dk_dH=-2.2*(1-1./H12)./H12.^3 -k./H12;

CD2=2./Ret .*k;
dCD2_dH=2./Ret .*dk_dH;
dCD2_dRet= -CD2./Ret;

end

