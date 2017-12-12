function [ CD2, dCD2_dH,dCD2_dRet ] = CD2lam(H12,Ret,H32,mod)
%CD2LAM     calculates 2*Cd and its the derivate in respect to the shape parameter H12 and Ret
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used


%nu=evalin('base','nu');
%H32=ET/U;
%Ret=T*U/nu;

if nargin==3
% "new" models from Drela 1989
    if H12 < 4
        CD2         =H32/Ret *  ( 0.207 + 0.00205*(4-H12)^5.5 ); 
        dCD2_dH     =H32/Ret *  ( - 0.011275*(4-H12)^4.5 );
    else    
        CD2         =H32/Ret *  ( 0.207 -0.0016*(H12-4)^2 /(1+0.002*(H12-4)^2) );
        dCD2_dH     =H32/Ret *  ( -0.0032*(H12-4)/(1+0.02*(H12-4)^2)...
                                    +0.0016*(H12-4)^2*(0.004*H12-0.016)/(1+0.002*(H12-4)^2)^2);
    end
    dCD2_dRet=-CD2/Ret;
else
% "old" models from Drela 1987
    if H12 < 4
        CD2         = H32/Ret *  ( 0.207 + 0.00205*(4-H12)^5.5 ); 
        dCD2_dH     = H32/Ret *  ( - 0.011275*(4-H12)^4.5 );
        % H32 eingesetzt
        %CD2=1/Ret *  ( 0.207 + 0.00205*(4-H12)^5.5 ) ...
        %           *  (1.515+0.076*(4-H12)^2 /H12) ;%d3 eingesetzt
    else    
        CD2         =H32/Ret *  ( 0.207 -0.003*(H12-4)^2 /(1+0.002*(H12-4)^2) );
        dCD2_dH     =H32/Ret *  ( -0.006*(H12-4)/(1+0.02*(H12-4)^2)...
                                  +0.003*(H12-4)^2*(0.004*H12-0.016)/(1+0.002*(H12-4)^2)^2);
        % H32 eingesetzt
        %CD2=1/Ret *( 0.207-0.003*(H12-4)^2 /(1+0.002*(H12-4)^2) )...
        %          *(1.515+0.04*(H12-4)^2/H12);
    end
    dCD2_dRet=-CD2/Ret;
end



end

