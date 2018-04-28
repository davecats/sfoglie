function [xseparation, xreattach] = FindSeparationLoc( x,H, iTran, HLmax,HTmax )
%FINDSEPERATIONLOC finds the seperation location with the criterion H12>Hkrit
%                  It has to be pointed out, that in the TURBULENT part the criterion 
%                  is NECESSARY but NOT SUFFICIENT. Inverse calculation is
%                  used while H12>HTmax but turbulent reattachement might
%                  occur before.


xseparation=NaN;
xreattach=NaN;

HL=H(1:iTran-1);
HT=H(iTran:end);

% look for seperated branch
isep =find(HL>HLmax ,1,'first');
if isempty(isep) % no laminar seperation
    % seperation might occur in turbulent region
   isep =find(HT>HTmax,1 ,'first') + iTran-1 ;
   w1= (H(isep)-HTmax)/(H(isep)-H(isep-1));
else % laminar seperation
   w1= (H(isep)-HLmax)/(H(isep)-H(isep-1)); 
end


if ~isempty(isep) % only if seperation occurs
    % get seperation location
    w2=1-w1;
    xseparation=x(isep-1)*w1+x(isep)*w2;
    
    % find possible reattachement
    
    if isep<iTran-1
        ire =find(HL(isep+1:end)<HLmax ,1 ,'first')+isep;
    else
        ire=[];
    end
    
    if ~isempty(ire) % laminar reattachement occurs
        w1= (H(ire)-HLmax)/(H(ire)-H(ire-1));
        w2=1-w1;
        xreattach=x(ire-1)*w1+x(ire)*w2;
    else
        
        % look for turbulent reattachement
        ire =find(H(isep+1:end)<HTmax,1 ,'first') +isep;
        if ~isempty(ire)
            w1= (H(ire)-HTmax)/(H(ire)-H(ire-1));
            w2=1-w1;
            xreattach=x(ire-1)*w1+x(ire)*w2;
        end
        
    end
    
    
    
end








end

