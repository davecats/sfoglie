function [ H32, dH32_dH] = H32lam(H12,mod)
%DELTA3LAM  calculates H32 and its the derivate in respect to the shape parameter H12
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used


if nargin==1
    mod=1;
end

H32=zeros(size(H12));
dH32_dH=zeros(size(H12));
if mod==1
    % "new" models from Drela 1989
    ind1=find(H12<=4.35);
    ind2=find(H12>4.35);
    H12TMP=H12;

    % H12<4.35      
    H12=H12TMP(ind1);
    K=(H12-4.35);
    K2=(H12+1);

    H32(ind1)      = ( 0.0111*K.^2./K2 - 0.0278*K.^3./K2 - 0.0002*(H12.*K).^2 +1.528);
    dH32_dH(ind1)  = ( 0.0222*K./K2 - 0.0111*K.^2./K2.^2 ...
                      -0.0834*K.^2./K2 + 0.0278*K.^3./K2.^2 ...
                      -0.0004*H12.*K.*(K + H12)  );


    % H12>4.35 
    H12=H12TMP(ind2);
    K=(H12-4.35);

    H32(ind2)     = ( 0.015*K.^2./H12 +1.528 );
    dH32_dH(ind2) = ( 0.03 *K./H12 - 0.015*K.^2./H12.^2 );

%---------------------------------------------------------------   
else
    % "old" models from Drela 1987
    ind1=find(H12<=4);
    ind2=find(H12>4);
    H12TMP=H12;

    % H12<4      
    H12=H12TMP(ind1);

    H32(ind1)      = ( 1.515 + 0.076*(4-H12).^2 ./H12 );
    dH32_dH(ind1)  = ( (-0.152*(4-H12)./H12)-0.076*(4-H12).^2./H12.^2 );


    % H12>4 
    H12=H12TMP(ind2);

    H32(ind2)     = (1.515 + 0.04*(H12-4).^2 ./H12 );
    dH32_dH(ind2) = ( (0.08*(H12-4)./H12) - 0.04*(H12-4).^2./H12.^2 );

end
    
end

