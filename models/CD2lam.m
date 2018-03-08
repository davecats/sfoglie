function [ CD2, dCD2_dH,dCD2_dRet ] = CD2lam(H12,Ret,mod)
%CD2LAM     calculates 2*Cd/HS and its the derivate in respect to the shape parameter H12 and Ret
%           mod: 1 - models from Drela 1989 (default)
%                2 - models from Drela 1987
%                3 - models from Eppler 1980

if nargin==2;
    mod=1;
end


CD2=zeros(size(H12));
dCD2_dH=zeros(size(H12));
if mod==1 % "new" models from Drela 1989     

    ind1=find(H12<=4);
    ind2=find(H12>4);
    RetTMP=Ret;
    H12TMP=H12;

    % H12 < 4      
    Ret=RetTMP(ind1);
    H12=H12TMP(ind1);

    CD2(ind1)     =1./Ret .*  ( 0.207 + 0.00205*(4-H12).^5.5 ); 
    dCD2_dH(ind1) =1./Ret .*  ( - 0.011275*(4-H12).^4.5 );

    % H12 > 4 
    Ret=RetTMP(ind2);
    H12=H12TMP(ind2);

    k1=(H12-4);
    k2=1+0.02*k1.^2;

    CD2(ind2)     =1./Ret .*  ( 0.207 -0.0016*k1.^2 ./k2 );
    dCD2_dH(ind2) =1./Ret .*  ( -0.0032*k1.*( 1./k2 - 0.02*k1.^2./k2.^2) );


    dCD2_dRet=-CD2./RetTMP;
%---------------------------------------------------------------
elseif mod==2  % "old" models from Drela 1987
    
    ind1=find(H12<=4);
    ind2=find(H12>4);
    RetTMP=Ret;
    H12TMP=H12;

    % H12 < 4      
    Ret=RetTMP(ind1);
    H12=H12TMP(ind1);

    CD2(ind1)     = 1./Ret .*  ( 0.207 + 0.00205*(4-H12).^5.5 ); 
    dCD2_dH(ind1) = 1./Ret .*  ( - 0.011275*(4-H12).^4.5 );

    % H12 > 4 
    Ret=RetTMP(ind2);
    H12=H12TMP(ind2);

    CD2(ind2)     =1./Ret .*  ( 0.207 -0.003*(H12-4).^2 ./(1+0.02*(H12-4).^2) );
    dCD2_dH(ind2) =1./Ret .*  ( -0.006*(H12-4)./(1+0.02*(H12-4).^2)...
                                  +0.003*(H12-4).^2.*(0.04*H12-0.16)./(1+0.02*(H12-4).^2).^2);

    dCD2_dRet=-CD2./RetTMP;

elseif mod==3 % models from Eppler 1980
    
    [H32, dH32_dH]=H32lam(H12);
    
    %CD2    = 2./Ret .*  (7.853976 - 10.260551*H32 + 3.418898* H32.^2 ); %
    CD2    = 2./Ret .* ( 7.853976./H32      - 10.260551 + 3.418898* H32 );
    dCD2_dH= 2./Ret .* (-7.853976./(H32.^2)             + 3.418898      ).*dH32_dH;
    
    dCD2_dRet=-CD2./Ret;
end
   


end

