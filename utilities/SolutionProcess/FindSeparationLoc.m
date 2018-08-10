function [xseparation, xreattach] = FindSeparationLoc( x,H, iTran, HLmax,HTmax,xtran )
%FINDSEPERATIONLOC finds the seperation location with the criterion H12>Hkrit
%                  It has to be pointed out, that in the TURBULENT part the criterion 
%                  is NECESSARY but NOT SUFFICIENT. Inverse calculation is
%                  used while H12>HTmax but turbulent reattachement might
%                  occur before.


xseparation=NaN;
xreattach=NaN;

HL=H(1:iTran-1);
HT=H(iTran:end);


% nodes with laminar separation
isep_L=find(HL>HLmax);
if ~isempty(isep_L) % separation occurs

    % weigthening factor
    w1= (H(isep_L(1))-HLmax)/(H(isep_L(1))-H(isep_L(1)-1)); 
    % secant method to find separation location
    xseparation= w1*x(isep_L(1)-1) + (1-w1)*x(isep_L(1));
    
    % look for laminar reattachement
    if isep_L(end)<iTran-1
        % weigthening factor
        w1= (HLmax-H(isep_L(end)+1))/(H(isep_L(end))-H(isep_L(end)+1));
        xreattach=w1*x(isep_L(end)) + (1-w1)*x(isep_L(end)+1);
    
    end
    
end


isep_T=find(HT>HTmax);
isep_T=isep_T+iTran-1;
if ~isempty(isep_T)
    % turbulent flow starts is not separated
    if isnan(xseparation) || (~isnan(xseparation) && ~isnan(xreattach)) 
        if isep_T(1)==iTran % if separation occurs right with transition
            xseparation=xtran;
        else
            % weigthening factor
            w1= (H(isep_T(1))-HTmax)/(H(isep_T(1))-H(isep_T(1)-1)); 
            % secant method to find separation location
            xseparation= w1*x(isep_T(1)-1) + (1-w1)*x(isep_T(1));
        end
    end
        % look for turbulent reattachement
    if isep_T(end)<length(H)
        % weigthening factor
        w1= (HTmax-H(isep_T(end)+1))/(H(isep_T(end))-H(isep_T(end)+1));
        xreattach=w1*x(isep_T(end)) + (1-w1)*x(isep_T(end)+1);

    end

end



end

