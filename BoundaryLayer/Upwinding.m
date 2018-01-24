function [ upw, dupw_dH1,dupw_dH2] = Upwinding(H1,H2,IsWake)
%UPWINDING Sets the upwind parameter vor shape parameter Equation
%          upw=0.5 -> centraldifferences
%          upw=1   -> forward differences


I=ones(length(H2),1);
if IsWake; k=1; else k=5; end
HC=k./H2.^2;
HL=log(abs((H2-I)./(H1-I)));
HLS=HL.^2; HLS(HLS>15)=15;
EHH=exp(-HLS.*HC);
upw=1-0.5*EHH;
    
dupw_dHL=EHH.*HL.*HC ;   
dupw_dHC=0.5*EHH.*HLS; 

dHL_dH1=-1./(H1-I);
dHL_dH2= 1./(H2-I);         

dHC_DH2=-2*HC./H2;

dupw_dH1=dupw_dHL.*dHL_dH1;
dupw_dH2=dupw_dHL.*dHL_dH2 + dupw_dHC.*dHC_DH2;

end

