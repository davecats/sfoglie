function [ Cf2, dCf2_dH, dCf2_dRet ] = CF2lam( H12, Ret, mod)
%CF2LAM     calculates Cf/2 and its the derivate in respect to the shape parameter H12 and Ret
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used
%   

% xtst=(0.1:0.1:5.5);
% xtst2=(5.5:0.1:10);
% I=ones(1,length(xtst));
% I2=ones(1,length(xtst2));
% a=  (  0.0727*(5.5*I-xtst).^3 ./(xtst+I) - 0.07*I );
% b=0.015*(I2-1./(xtst2-4.5*I2)).^2 -0.07*I2;
% figure
% hold on
% plot(xtst,a);
% plot(xtst2,b);




% if vector input
%---------------------------------------------------------------
if isscalar(H12)==false
    Cf2=zeros(size(H12));
    dCf2_dH=zeros(size(H12));
    if nargin==2
        % "new" models from Drela 1989
        %ind1=find(H12<5.5 & Ret~=0);
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
    else
        % "old" models from Drela 1987
        ind1=find(H12<7.4);
        ind2=find(H12>=7.4);
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<7.4      
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
        
        Cf2(ind1)     = 1./Ret.* ( -0.067 + 0.01977*(7.4-H12).^2 ./(H12-1) ); 
        dCf2_dH(ind1) = 1./Ret.* ( -0.03954*(7.4-H12)./(H12-1) ...
                                   -0.01977*(7.4-H12).^2./(H12-1).^2 ); 


        % H12>7.4 
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        Cf2(ind2)     = 1./Ret.* ( -0.067 + 0.022*(1-1.4./(H12-6)).^2 );
        dCf2_dH(ind2) = 1./Ret.* ( 0.0616* (1-1.4./(H12-6))./(H12-6).^2  );
            
       dCf2_dRet=-Cf2./RetTMP;      
    end
    
     % singularity
     %dCf2_dRet(RetTMP==0)=0;
    
     
     
% if scalar input
%---------------------------------------------------------------
else
        %singularity
    if Ret==0
        Cf2=0;
        dCf2_dH=0;
        dCf2_dRet=0;
    elseif nargin==2
    % "new" models from Drela 1989
         if  H12< 5.5
            Cf2     =  1/Ret* (  0.0727*(5.5-H12)^3 /(H12+1) - 0.07 ); 
            dCf2_dH =  1/Ret* ( -0.2181*(5.5-H12)^2 /(H12+1) - 0.0727*(5.5-H12)^3 /(H12+1)^2 );
         else
            Cf2     = 1/Ret* ( 0.015*(1-1/(H12-4.5))^2 -0.07 ); 
            dCf2_dH = 1/Ret*  ( 0.03* (1-1.4/(H12-6))/(H12-6)^2  );
         end
         dCf2_dRet= -Cf2/Ret;
         
         
         % factor for cf/2
         Cf2= Cf2/2;
         dCf2_dH=dCf2_dH/2;
         dCf2_dRet=dCf2_dRet/2;
         
         
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

