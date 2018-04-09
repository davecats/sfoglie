function [ Cf2,dCf2_dH,dCf2_dRet ] = CF2turb( H,Ret ,onlyTurb,mod)
%CF2TURB  calculates Cf/2 and its the derivate in respect to the shape
%         parameter H12 and Ret in turbulent case
%         returns the maximum value of laminar and turbulent Cf. If onlyTurb=true do not check for bigger laminar Cf  


if nargin<4
    mod=1;
end

if mod==1 || mod==2 %  models from Drela 1987/89
    
    e1=-1.33*H; 
    e1(e1<-20)=-20;
    e2=-1.74 - 0.31*H;

    %tst= sin(0.0048*Ret).*max(370-Ret,0)./(370-Ret) + (sin(0.0048*370)*(exp(-0.016.*(Ret-370)))).*max(Ret-370,0)./(Ret-370);
    
    
    logRet= log(Ret); 
    logRet(logRet<3)=3;

    T1=0.3*exp(e1).*(logRet/2.3025851).^e2;
    TH=tanh(4 - H/0.875);
    T2=0.00011*(TH-1);

    Cf2= T1 + T2 ;%+ 0.0016*tst;
    dCf2_dH   = ( -1.33 -0.31*log(logRet/2.3025851) ) .*T1 ...
                -0.00011*(1-TH.^2)/0.875;


    dCf2_dRet = e2.*T1./(logRet.*Ret);

    % factor for cf/2
    Cf2=0.5*Cf2;
    dCf2_dH   =0.5*dCf2_dH ;
    dCf2_dRet =0.5*dCf2_dRet ; 
    
elseif mod==3  % models from Eppler 1980
    
    e1=-1.26*H; 
    %e1(e1<-20)=-20;
    
    Cf2      = 0.045716*( (H-1).*Ret ).^(-0.232).*exp( e1 );
    dCf2_dH  =-0.010606112*(H-1).^(-1.232).*Ret.^(-0.232).*exp( e1 ) - 0.05760216*( (H-1).*Ret ).^(-0.232).*exp( e1 );
    dCf2_dRet=-0.010606112*(H-1).^(-0.232).*Ret.^(-1.232).*exp( e1 );
end

if nargin==2 || onlyTurb==false
    % if laminar Cf is bigger -> use laminar value
    [CFL,~,~]=CF2lam( H,Ret,mod );
    ind= find(CFL>Cf2);
    if ~isempty(ind)
        [CFL,dCFL_dH,dCFL_dRet]=CF2lam( H(ind),Ret(ind),mod );
        Cf2(ind)=CFL;
        dCf2_dH(ind)  = dCFL_dH;
        dCf2_dRet(ind)= dCFL_dRet;
        %disp('CF set');
    end
end



end

