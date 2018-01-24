function [ out ] = AngleIntervall( in )
%ANGLEINTERVALL     transforms an angle between -pi and pi to one between 0 and 2pi
%                   vice versa

% angle between -pi and pi
if isempty(find(in<0,1))== false
    out=in;
    out(in<0)=2*pi+in(in<0);
% angle between 0 and 2pi
else
    out=in;
    out(in>pi)=in(in>pi)-2*pi;
end


end

