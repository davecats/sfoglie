function [ Nle,sLE,LE1,LE2 ] = getStagnationPoint( gam, s )
%GETSTAGNATIONPOINT find stagnation point location with secant approximation


Nle=find(gam<0,1); % first node of pressure side
 
 L=s(Nle)-s(Nle-1);

 if gam(Nle-1) < -gam(Nle)
     sLE= s(Nle-1) - L*gam(Nle-1)/(gam(Nle)-gam(Nle-1));
 else
     sLE= s(Nle) - L*gam(Nle)/(gam(Nle)-gam(Nle-1));
 end
 
 % if stagnation point is close to a node -> make shure it falls in intervall
 if sLE < s(Nle-1); sLE = s(Nle-1) +1e-7; end
 if sLE > s(Nle  ); sLE = s(Nle  ) -1e-7; end
 
 LE1= sLE - s(Nle-1); % distance from last suction side point
 LE2 =L-LE1;% distance from first pressure side point

end

