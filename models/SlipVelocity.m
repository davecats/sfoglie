function [Us,dUs_dH12,dUs_dHS] = SlipVelocity(H12,HS)
%SLIPVELOCITY calculates the velocity at the edge of the inner layer Us


Us=0.5*HS.*( -1/3 + 1./(0.75*H12) );

dUs_dHS=  Us./HS;
dUs_dH12 = -HS./(1.5*H12.^2);

end

