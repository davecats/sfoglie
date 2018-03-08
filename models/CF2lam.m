function [ Cf2, dCf2_dH, dCf2_dRet ] = CF2lam( H12, Ret, mod)
%CF2LAM     calculates Cf/2 and its the derivate in respect to the shape parameter H12 and Ret
%           mod: 1 - models from Drela 1989 (default)
%                2 - models from Drela 1987
%                3 - models from Eppler 1980

if nargin==2;
    mod=1;
end

Cf2=zeros(size(H12));
dCf2_dH=zeros(size(H12));
if mod==1 % "new" models from Drela 1989
    ind1=find(H12<5.5 );
    ind2=find(H12>=5.5);
    RetTMP=Ret;
    H12TMP=H12;

    % H12<5.5     
    Ret=RetTMP(ind1);
    H12=H12TMP(ind1);

    Cf2(ind1)     =  1./Ret .* (  0.0727*(5.5-H12).^3 ./(H12+1) - 0.07 ); 
    dCf2_dH(ind1) =  1./Ret .* ( -0.2181*(5.5-H12).^2 ./(H12+1) ...
                                 - 0.0727*(5.5-H12).^3 ./(H12+1).^2 );
    % H12>5.5
    Ret=RetTMP(ind2);
    H12=H12TMP(ind2);

    Cf2(ind2)      = 1./Ret.* ( 0.015*(1-1./(H12-4.5)).^2 -0.07 ); 
    dCf2_dH(ind2)  = 1./Ret.* ( 0.03* (1-1./(H12-4.5))./(H12-4.5).^2  );


    dCf2_dRet=-Cf2./RetTMP;

    % factor for cf/2
    dCf2_dRet=dCf2_dRet/2;
    dCf2_dH=dCf2_dH/2;
    Cf2=Cf2/2;
    %---------------------------------------------------------------
elseif mod==2  % "old" models from Drela 1987
    ind1=find(H12<7.4);
    ind2=find(H12>=7.4);
    RetTMP=Ret;
    H12TMP=H12;

    % H12 < 7.4      
    Ret=RetTMP(ind1);
    H12=H12TMP(ind1);

    Cf2(ind1)     = 1./Ret.* ( -0.067 + 0.01977*(7.4-H12).^2 ./(H12-1) ); 
    dCf2_dH(ind1) = 1./Ret.* ( -0.03954*(7.4-H12)./(H12-1) ...
                               -0.01977*(7.4-H12).^2./(H12-1).^2 ); 


    % H12 > 7.4 
    Ret=RetTMP(ind2);
    H12=H12TMP(ind2);

    Cf2(ind2)     = 1./Ret.* ( -0.067 + 0.022*(1-1.4./(H12-6)).^2 );
    dCf2_dH(ind2) = 1./Ret.* ( 0.0616* (1-1.4./(H12-6))./(H12-6).^2  );

    dCf2_dRet=-Cf2./RetTMP;    
elseif mod==3 % models from Eppler 1980
    
    [H32, dH32_dH]=H32lam(H12);
    
    ind1=find(H32>=1.51509 & H32<1.57258);
    ind2=find(H32>=1.57258);

    % 1.51509 < H32 < 1.57258
    
    Cf2(ind1)    = 1./Ret(ind1).*( 2.512589 - 1.686095*H12(ind1)  + 0.391541*H12(ind1).^2 - 0.031729*H12(ind1).^3 );
    dCf2_dH(ind1)= 1./Ret(ind1).*(          - 1.686095            + 0.783082*H12(ind1)    - 0.095187*H12(ind1).^2 );
    
    
    % H32 > 1.57258
    Cf2(ind2)    = 1./Ret(ind2).*( 1.372391 - 4.226253*H32(ind2) + 2.221687*H32(ind2).^2 );
    dCf2_dH(ind2)= 1./Ret(ind2).*(          - 4.226253           + 4.443374*H32(ind2)    ).*dH32_dH(ind2);
    
    dCf2_dRet=-Cf2./Ret;  
end
   
   
end

