function Rel = GetRelaxationFactor( delta, phi )
%GETRELAXATIONFACTOR calculate factor for underrelaxation if needed

Rel=ones(size(phi));

Rel(delta./phi>1.5)=1.5*phi(delta./phi>1.5)./delta(delta./phi>1.5);
Rel(delta./phi<-0.5)=-0.5*phi(delta./phi<-0.5)./delta(delta./phi<-0.5);


end

