function [ ET, dET_dH] = delta3lam(T,H12,mod)
%DELTA3LAM  calculates energy thickness and its the derivate in respect to the shape parameter H12
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used


if isscalar(H12)==false
    n=length(H12);
    ET=zeros(n,1);
    dET_dH=zeros(n,1);
    if nargin==2
        % "new" models from Drela 1989
        ind1=find(H12<=4.35);
        ind2=find(H12>4.35);
        TTMP=T;
        H12TMP=H12;
        
        % H12<4      
        I=ones(length(ind1),1);
        T=TTMP(ind1);
        H12=H12TMP(ind1);
        
        ET(ind1)      = T.*( 0.0111*K.^2./K2 - 0.0278*K.^3./K2-0.0002*(H12.*K).^2 +1.528*I);
        dET_dH(ind1)  = T.*( 0.0222*K./K2-0.0111*K.^2./K2.^2 ...
                            -0.0834*K.^2./K2 + 0.0278*K.^3./K2.^2 ...
                            -0.0004*H12.*K.*(K+H12)  );
        
        
        % H12>4 
        I=ones(length(ind2),1);
        T=TTMP(ind2);
        H12=H12TMP(ind2);
        
        ET(ind2)     = T.*( 0.015*K.^2./H12 +1.528 );
        dET_dH(ind2) = T.*( 0.03*K./H12 - 0.015*K.^2./H12.^2 );
       
    else
        % "old" models from Drela 1987
        ind1=find(H12<=4);
        ind2=find(H12>4);
        TTMP=T;
        H12TMP=H12;
        
        % H12<4      
        I=ones(length(ind1),1);
        T=TTMP(ind1);
        H12=H12TMP(ind1);
        
        
        ET(ind1)      = T.*( 1.515*I + 0.076*(4*I-H12).^2 ./H12 );
        dET_dH(ind1)  = T.*( (-0.152*(4*I-H12)./H12)-0.076*(4*I-H12).^2./H12.^2 );



        % H12>4 
        I=ones(length(ind2),1);
        T=TTMP(ind2);
        H12=H12TMP(ind2);
        
        ET(ind2)     = T.*(1.515*I + 0.04*(H12-4*I).^2 ./H12 );
        dET_dH(ind2) = T.*( (0.08*(H12-4*I)./H12) - 0.04*(H12-4*I).^2./H12.^2 );
        
    end
    

else
    if nargin==2
    % "new" models from Drela 1989
         K=(H12-4.35); 
         if H12 < 4.35
            ET      =T* ( 0.0111*K^2/(H12+1) - 0.0278*K^3/(H12+1)-0.0002*(H12*K)^2 +1.528);
            dET_dH  =T* ( 0.0222*K/(H12+1)-0.0111*K^2/(H12+1)^2 ...
                         -0.0834*K^2/(H12+1) + 0.0278*K^3/(H12+1)^2 ...
                         -0.0004*H12*K*(K+H12)  );
         else    
            ET     =T* ( 0.015*K^2/H12 +1.528 );
            dET_dH =T* ( 0.03*K/H12 - 0.015*K^2/H12^2 );
         end
    else
        % "old" models from Drela 1987
         if H12 < 4
            ET      =T* ( 1.515 + 0.076*(4-H12)^2 /H12 );
            dET_dH  =T* ( (-0.152*(4-H12)/H12)-0.076*(4-H12)^2/H12^2 );
         else    
            ET     =T* (1.515 + 0.04*(H12-4)^2 /H12 );
            dET_dH =T* ( (0.08*(H12-4)/H12) - 0.04*(H12-4)^2/H12^2 );
         end
    end
end

end

