function [ CtEQ,dCt_dH,dCt_dRet,dCt_dHS, dCt_dUs] = CtEQ( H,Ret,HS,Us,IsWake)
%CTEQ Calculates equilibrium shear stress coefficient

% CtEQ
%--------------------------------------------------------------
dim=size(H);

dk_dH=ones(dim);
if IsWake
    k= H - 1 ;
    dk_dRet=zeros(dim);

else
    k= H - 1 - 18./Ret;
    dk_dRet=18./Ret.^2;
    
    ind2=find(k<0.01);
    %if ~isempty(ind2); disp('in CtEQ: k<0,01') ;end
    k(ind2)=0.01;
    dk_dRet(ind2)=0;
    dk_dH(ind2)=0;
end

H(H<1.00005)=1.00005;
Us(Us>0.99999)=0.99999;
k2=0.014851118*HS.*(H-1).*k.^2 ./ ((1-Us).*H.^3);

dk2_dUs=k2./(1-Us);
dk2_dHS=k2./HS;
dk2_dk=2*k2./k;

CtEQ= sqrt( k2 );
dCt_dk2= 1./(2*CtEQ);

dCt_dUs=dCt_dk2.*dk2_dUs;

dCt_dHS=dCt_dk2.*dk2_dHS ;%+ dCt_dUs.*dUs_dHS;

dCt_dRet=dCt_dk2 .*dk2_dk .*dk_dRet ;

dCt_dH=dCt_dk2 .* (k2./(H-1) - 3*k2./H + dk2_dk.*dk_dH ) ;%+ dCt_dUs.*dUs_dH ;


end

