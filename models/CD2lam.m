function [ CD2, dCD2_dH,dCD2_dRet ] = CD2lam(H12,Ret,H32,mod)
%CD2LAM     calculates 2*Cd and its the derivate in respect to the shape parameter H12 and Ret
%           default: takes the "new" models
%           if mod is committed, the models of the paper are used

% if vector input
%---------------------------------------------------------------
if isscalar(H12)==false
    n=length(H12);
    CD2=zeros(n,1);
    dCD2_dH=zeros(n,1);
    if nargin==3
        % "new" models from Drela 1989
        ind1=find(H12<=4);
        ind2=find(H12>4);
        H32TMP=H32;
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<4      
        I=ones(length(ind1),1);
        H32=H32TMP(ind1);
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
        
        CD2(ind1)     =H32./Ret .*  ( 0.207*I + 0.00205*(4*I-H12).^5.5 ); 
        dCD2_dH(ind1) =H32./Ret .*  ( - 0.011275*(4*I-H12).^4.5 );
        
        % H12>4 
        I=ones(length(ind2),1);
        H32=H32TMP(ind2);
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        
        CD2(ind2)     =H32./Ret .*  ( 0.207*I -0.0016*(H12-4*I).^2 ./(I+0.002*(H12-4*I).^2) );
        dCD2_dH(ind2) =H32./Ret .*  ( -0.0032*(H12-4*I)./(I+0.02*(H12-4*I).^2)...
                                      +0.0016*(H12-4*I).^2.*(0.004*H12-0.016*I)./(I+0.002*(H12-4*I).^2).^2);
        
        
        
        dCD2_dRet=-CD2./RetTMP;
    %---------------------------------------------------------------
    else
        % "old" models from Drela 1987
        ind1=find(H12<=4);
        ind2=find(H12>4);
        H32TMP=H32;
        RetTMP=Ret;
        H12TMP=H12;
        
        % H12<4      
        I=ones(length(ind1),1);
        H32=H32TMP(ind1);
        Ret=RetTMP(ind1);
        H12=H12TMP(ind1);
    
        CD2(ind1)     = H32./Ret .*  ( 0.207*I + 0.00205*(4*I-H12).^5.5 ); 
        dCD2_dH(ind1) = H32./Ret .*  ( - 0.011275*(4*I-H12).^4.5 );

        % H12>4 
        I=ones(length(ind2),1);
        H32=H32TMP(ind2);
        Ret=RetTMP(ind2);
        H12=H12TMP(ind2);
        
        CD2(ind2)     =H32./Ret .*  ( 0.207*I -0.003*(H12-4*I).^2 ./(I+0.002*(H12-4*I).^2) );
        dCD2_dH(ind2) =H32./Ret .*  ( -0.006*(H12-4*I)./(I+0.02*(H12-4*I).^2)...
                                      +0.003*(H12-4*I).^2.*(0.004*H12-0.016*I)./(I+0.002*(H12-4*I).^2).^2);

        dCD2_dRet=-CD2./RetTMP;
    end
    
% if scalar input
%---------------------------------------------------------------
else
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




end

