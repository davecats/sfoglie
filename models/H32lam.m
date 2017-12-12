function [ H32, dH32_dH] = H32lam(H12,mod)
%DELTA3LAM  calculates H32 and its the derivate in respect to the shape parameter H12
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used

if nargin==1
% "new" models from Drela 1989
     K=(H12-4.35); 
     if H12 < 4.35
        H32      = ( 0.0111*K^2/(H12+1) - 0.0278*K^3/(H12+1)-0.0002*(H12*K)^2 +1.528);
        dH32_dH  = ( 0.0222*K/(H12+1)-0.0111*K^2/(H12+1)^2 ...
                     -0.0834*K^2/(H12+1) + 0.0278*K^3/(H12+1)^2 ...
                     -0.0004*H12*K*(K+H12)  );
     else    
        H32     = ( 0.015*K^2/H12 +1.528 );
        dH32_dH = ( 0.03*K/H12 - 0.015*K^2/H12^2 );
     end
else
    % "old" models from Drela 1987
     if H12 < 4
        H32      = ( 1.515 + 0.076*(4-H12)^2 /H12 );
        dH32_dH  = ( (-0.152*(4-H12)/H12)-0.076*(4-H12)^2/H12^2 );
     else    
        H32     = (1.515 + 0.04*(H12-4)^2 /H12 );
        dH32_dH = ( (0.08*(H12-4)/H12) - 0.04*(H12-4)^2/H12^2 );
     end
end
 
end

