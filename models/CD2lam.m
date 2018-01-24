function [ CD2, dCD2_dH,dCD2_dRet ] = CD2lam(H12,Ret,mod)
%CD2LAM     calculates 2*Cd/HS and its the derivate in respect to the shape parameter H12 and Ret
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used


% xtst=(0.1:0.1:4);
% xtst2=(4:0.1:10);
% I=ones(1,length(xtst));
% I2=ones(1,length(xtst2));
% a=  ( 0.207*I + 0.00205*(4*I-xtst).^5.5 );
% b= ( 0.207*I2 -0.0016*(xtst2-4*I2).^2 ./(I2+0.02*(xtst2-4*I2).^2) );
% figure
% hold on
% plot(xtst,a);
% plot(xtst2,b);




% if vector input
%---------------------------------------------------------------
if isscalar(H12)==false
    CD2=zeros(size(H12));
    dCD2_dH=zeros(size(H12));
    if nargin==2
        % "new" models from Drela 1989        
        
        %ind1=find(H12<=4 & Ret~=0);
        ind1=find(H12<=4);
        ind2=find(H12>4);
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<4      
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
        
        CD2(ind1)     =1./Ret .*  ( 0.207 + 0.00205*(4-H12).^5.5 ); 
        dCD2_dH(ind1) =1./Ret .*  ( - 0.011275*(4-H12).^4.5 );
        
        % H12>4 
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        k1=(H12-4);
        k2=1+0.02*k1.^2;
        
        CD2(ind2)     =1./Ret .*  ( 0.207 -0.0016*k1.^2 ./k2 );
        dCD2_dH(ind2) =1./Ret .*  ( -0.0032*k1.*( 1./k2 - 0.04*k1.^2./k2.^2) );
        

        dCD2_dRet=-CD2./RetTMP;
    %---------------------------------------------------------------
    else
        % "old" models from Drela 1987
        ind1=find(H12<=4);
        ind2=find(H12>4);
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<4      
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
    
        CD2(ind1)     = 1./Ret .*  ( 0.207 + 0.00205*(4-H12).^5.5 ); 
        dCD2_dH(ind1) = 1./Ret .*  ( - 0.011275*(4-H12).^4.5 );

        % H12>4 
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        CD2(ind2)     =1./Ret .*  ( 0.207 -0.003*(H12-4).^2 ./(1+0.02*(H12-4).^2) );
        dCD2_dH(ind2) =1./Ret .*  ( -0.006*(H12-4)./(1+0.02*(H12-4).^2)...
                                      +0.003*(H12-4).^2.*(0.04*H12-0.16)./(1+0.02*(H12-4).^2).^2);

        dCD2_dRet=-CD2./RetTMP;
    end
    
%     % singularity
%     dCD2_dRet(RetTMP==0)=0;

% if scalar input
%---------------------------------------------------------------
else

    %singularity
    if Ret==0
        CD2=0;
        dCD2_dH=0;
        dCD2_dRet=0;
    elseif nargin==2 
    % "new" models from Drela 1989
        if H12 < 4
            CD2         =1/Ret *  ( 0.207 + 0.00205*(4-H12)^5.5 ); 
            dCD2_dH     =1/Ret *  ( - 0.011275*(4-H12)^4.5 );
        else    
            CD2         =1/Ret *  ( 0.207 -0.0016*(H12-4)^2 /(1+0.02*(H12-4)^2) );
            dCD2_dH     =1/Ret *  ( -0.0032*(H12-4)/(1+0.02*(H12-4)^2)...
                                        +0.0016*(H12-4)^2*(0.04*H12-0.16)/(1+0.02*(H12-4)^2)^2);
        end
        dCD2_dRet=-CD2/Ret;
    else
    % "old" models from Drela 1987
        if H12 < 4
            CD2         = 1/Ret *  ( 0.207 + 0.00205*(4-H12)^5.5 ); 
            dCD2_dH     = 1/Ret *  ( - 0.011275*(4-H12)^4.5 );
        else    
            CD2         =1/Ret *  ( 0.207 -0.003*(H12-4)^2 /(1+0.02*(H12-4)^2) );
            dCD2_dH     =1/Ret *  ( -0.006*(H12-4)/(1+0.02*(H12-4)^2)...
                                      +0.003*(H12-4)^2*(0.04*H12-0.16)/(1+0.02*(H12-4)^2)^2);

        end
        dCD2_dRet=-CD2/Ret;
    end
end




end

