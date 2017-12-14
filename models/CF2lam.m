function [ Cf2, dCf2_dH, dCf2_dRet ] = CF2lam( H12, Ret, mod)
%CF2LAM     calculates Cf/2 and its the derivate in respect to the shape parameter H12 and Ret
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used
%   

% if vector input
%---------------------------------------------------------------
if isscalar(H12)==false
    n=length(H12);
    Cf2=zeros(n,1);
    dCf2_dH=zeros(n,1);
    dCf2_dRet=zeros(n,1);
    if nargin==2
        % "new" models from Drela 1989
        ind1=find(H12<=5.5);
        ind2=find(H12>5.5);
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<5.5     
        I=ones(length(ind1),1);
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
        
        Cf2(ind1)     =  1./Ret .* (  0.0727*(5.5*I-H12).^3 ./(H12+I) - 0.07*I ); 
        dCf2_dH(ind1) =  1./Ret .* ( -0.2181*(5.5*I-H12).^2 ./(H12+I) ...
                                     - 0.0727*(5.5*I-H12).^3 ./(H12+I).^2 );
        
        % H12>5.5
        I=ones(length(ind2),1);
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        Cf2(ind2)      = 1./Ret.* ( 0.015*(I-1./(H12-4.5*I)).^2 -0.07*I ); 
        dCf2_dH(ind2)  = 1./Ret.*  ( 0.03* (I-1.4./(H12-6*I))./(H12-6*I).^2  );
       
        
        dCf2_dRet=-Cf2./RetTMP;
    %---------------------------------------------------------------
    else
        % "old" models from Drela 1987
        ind1=find(H12<=7.4);
        ind2=find(H12>7.4);
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<7.4      
        I=ones(length(ind1),1);
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
        
        Cf2(ind1)     = 1./Ret.* ( -0.067*I + 0.01977*(7.4*I-H12).^2 ./(H12-I) ); 
        dCf2_dH(ind1) = 1./Ret.* ( -0.03954*(7.4*I-H12)./(H12-I) ...
                                   -0.01977*(7.4*I-H12).^2./(H12-I).^2 ); 


        % H12>7.4 
        I=ones(length(ind2),1);
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        Cf2(ind2)     = 1./Ret.* ( -0.067*I + 0.022*(I-1.4./(H12-6*I)).^2 );
        dCf2_dH(ind2) = 1./Ret.* ( 0.0616* (I-1.4./(H12-6*I))./(H12-6*I).^2  );
            
       dCf2_dRet=-Cf2./RetTMP;      
    end
    
% if scalar input
%---------------------------------------------------------------
else
    %nu=evalin('base','nu');
    %Ret=(U*T)/nu;

    if nargin==2
    % "new" models from Drela 1989
         if  H12< 5.5
            Cf2     =  1/Ret* (  0.0727*(5.5-H12)^3 /(H12+1) - 0.07 ); 
            dCf2_dH =  1/Ret* ( -0.2181*(5.5-H12)^2 /(H12+1) - 0.0727*(5.5-H12)^3 /(H12+1)^2 );
         else
            Cf2     = 1/Ret* ( 0.015*(1-1/(H12-4.5))^2 -0.07 ); 
            dCf2_dH = 1/Ret*  ( 0.03* (1-1.4/(H12-6))/(H12-6)^2  );
         end
         dCf2_dRet= -Cf2/Ret;
    else
    % "old" models from Drela 1987
         if  H12< 7.4
            Cf2     = 1/Ret* ( -0.067 + 0.01977*(7.4-H12)^2 /(H12-1) ); 
            dCf2_dH = 1/Ret* ( -0.03954*(7.4-H12) /(H12-1) - 0.01977*(7.4-H12)^2/(H12-1)^2 ); 
         else
            Cf2     = 1/Ret* ( -0.067 + 0.022*(1-1.4/(H12-6))^2 );
            dCf2_dH = 1/Ret* ( 0.0616* (1-1.4/(H12-6))/(H12-6)^2  );
         end
         dCf2_dRet= -Cf2/Ret;
    end
end

end

