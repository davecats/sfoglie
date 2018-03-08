function [ CD2, dCD_dH,dCD_dRet, dCD_dHS,dCD_dUs, dCD_dCt] = CD2turb( H,Ret,HS,Us,Ct,IsWake)
%CD2TURB  calculates 2CD/HS aswell as sqrt(CtauEQ) and its the derivate in respect to the shape
%         parameter H12, Ret, H32,Cf, and Ct in turbulent case


% only wall contribution
dim=size(H);

%dCD_dCt=zeros(dim);
dCD_dH =zeros(dim);

if IsWake==false
    [ Cf,dCf_dH,dCf_dRet ]=CF2turb( H,Ret,true );
    CD2=2* Cf.*Us./HS;

    dCD_dUs=  CD2./Us;
    dCD_dHS= -CD2./HS;
    dCD_dCf=  CD2./Cf;
    
    dCD_dH  =dCD_dCf.*dCf_dH;
    dCD_dRet=dCD_dCf.*dCf_dRet;
    
    %correction factor for small H
    logR=log(Ret);
    Hmin=1 + 2.1 ./logR;
    dHmin_Ret=-2.1 ./ (logR.^2.*Ret);
    F=(H-1)./(Hmin-1);
    dF_dH= 1./(Hmin-1);
    dF_dRet= -F./(Hmin-1) .*dHmin_Ret;
    TF=tanh(F);
    FAC=0.5 + 0.5*TF;
    dFAC_dF  = 0.5*(1-TF.^2);
    dFAC_dH  =dFAC_dF.*dF_dH;
    dFAC_dRet=dFAC_dF.*dF_dRet;
    
    dCD_dH =FAC.*dCD_dH + CD2.*dFAC_dH ;
    dCD_dRet=FAC.*dCD_dRet + dFAC_dRet.*CD2;
    dCD_dHS=FAC.*dCD_dHS;
    dCD_dUs=FAC.*dCD_dUs;
    dCD_dCf=FAC.*dCD_dCf;
    %dCD_dCt=FAC.*dCD_dCt;
    
    CD2=FAC.*CD2;
else % wake -> falls away because Cf=0 on wake
    CD2=zeros(dim);
    dCD_dUs= zeros(dim);
    dCD_dHS= zeros(dim);
    dCD_dCf= zeros(dim);
    dCD_dRet=zeros(dim);
end


% add outer layer contribution
add=2*Ct.^2 .*(0.995-Us)./HS;
dadd_dUs= -2*Ct.^2 ./HS;
dadd_dHS= -add./HS;
dadd_dCt=  2*add./Ct;

CD2= CD2 + add;
dCD_dHS= dCD_dHS + dadd_dHS;
dCD_dUs= dCD_dUs + dadd_dUs;
%dCD_dCt= dCD_dCt + dadd_dCt;
dCD_dCt= dadd_dCt;

% add laminar stress contribution
add=0.3*(0.995-Us).^2./(HS.*Ret);
dadd_dUs = -2*add./(0.995-Us);
dadd_dHS = - add./HS;
dadd_dRet= - add./Ret;

CD2= CD2 + add;
dCD_dHS = dCD_dHS + dadd_dHS ;
dCD_dUs = dCD_dUs + dadd_dUs ;
dCD_dRet= dCD_dRet+ dadd_dRet;

% if laminar model gives biger values take those
if IsWake==false
   [CDL,~,~]=CD2lam( H,Ret );
    ind= find(CDL>CD2);
    if ~isempty(ind)
        [CDL,dCDL_dH,dCDL_dRet]=CD2lam( H(ind),Ret(ind) );
        CD2(ind)=CDL;
        dCD_dH(ind)  = dCDL_dH;
        dCD_dRet(ind)= dCDL_dRet;
        dCD_dHS(ind)=0;
        dCD_dUs(ind)=0;
        dCD_dCt(ind)=0;
        %dCD_dCf(ind)=0;
        %disp('CD set lam');
    end
else % Wake treatment
    [CDL,~,~]=CD2Lwake( H,Ret );
    ind= find(CDL>CD2);
    if ~isempty(ind)
        [CDL,dCDL_dH,dCDL_dRet]=CF2lam( H(ind),Ret(ind) );
        CD2(ind)=CDL;
        dCD_dH(ind)  = dCDL_dH;
        dCD_dRet(ind)= dCDL_dRet;
        dCD_dHS(ind)=0;
        dCD_dUs(ind)=0;
        dCD_dCt(ind)=0;
        %dCD_dCf(ind)=0;
        %disp('Wake: CD set lam');
   end
end

% Eppler 1980 -> without shear stress Equation
% CD2      = 0.01*( (H-1).*Ret ).^(-1/6);
% dCD2_dH  = -0.01/6*(H-1).^(-7/6).*Ret.^(-1/6);
% dCD2_dRet= -0.01/6*(H-1).^(-1/6).*Ret.^(-7/6);
end

