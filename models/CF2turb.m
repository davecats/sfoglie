function [ Cf2,dCf2_dH,dCf2_dRet ] = CF2turb( H,Ret )
%CF2TURB  calculates Cf/2 and its the derivate in respect to the shape
%         parameter H12 and Ret in turbulent case


e1=-1.33*H; 
e1(e1<-20)=-20;
e2=-1.74 - 0.31*H;
logRet= log(Ret); 
logRet(logRet<3)=3;

T1=0.3*exp(e1).*(logRet/2.3026).^e2;
TH=tanh(4 - H/0.875);
T2=0.00011*(TH-1);

Cf2= T1 + T2;
dCf2_dH   = ( -1.33*T1 -0.31*log(logRet/2.3026) ) .*T1 ...
            -0.00011*(1-TH.^2)/0.875;
      
        
dCf2_dRet = e2.*T1./(logRet.*Ret);

% factor for cf/2
Cf2=0.5*Cf2;
dCf2_dH =0.5*dCf2_dH ;
dCf2_dRet =0.5*dCf2_dRet ;


% if laminar Cf is bigger -> use laminar value
[CFL,~,~]=CF2lam( H,Ret );
ind= find(CFL>Cf2);
if ~isempty(ind)
    [CFL,dCFL_dH,dCFL_dRet]=CF2lam( H(ind),Ret(ind) );
    Cf2(ind)=CFL;
    dCf2_dH(ind)  = dCFL_dH;
    dCf2_dRet(ind)= dCFL_dRet;
    %disp('CF set');
end




end

