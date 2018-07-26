function [f,fR,fp] = StressVector(p,tau,n,Vb )
%STRESSVECTOR Calculates the stress vector for each airfoil node / panel using lemma of Cauchy.
%             input  
%             p:   pressure (colum vector)
%             tau: shear stresses  (colum vector)
%             n:   normal vector components (Matrix with 1 and 2 components)
%             Vb:  vector of wallnormal velocity (zero if not comitted)  
%             -> for dimensionless values the input musst already be normalized 
%           
%             output
%             f =[fx ,fy ] stress vector
%             fR=[fRx,fRy] contribution of shear stresses
%             fp=[fRx,fRy] contribution of pressure stresses



m=length(p);

% transpose n to have right dimensions if needed
if length(n(:,1))~=m
   if length(n(1,:))==m
        n=n'; 
   else
       err('Error in StrainVector: dimensions of input do not fit')
   end
end

% add force of blowing / suction in case Vb is committed
if nargin==4
    p = p + sign(Vb(1:m)).*Vb(1:m).^2;   
end

% find first point of
Nle=find(n(:,2)>0,1);

% switch sign for suction side -> force is in opposit direction of panel tangent vector
tau(1:Nle-1)=-tau(1:Nle-1);

%% Lemma von Cauchy

% in x-direction
fRx=  tau.*n(:,2); % shear force
fpx=  p  .*n(:,1); % pressure force

% in y-direction
fRy= -tau.*n(:,1); % shear force
fpy=  p  .*n(:,2); % pressure force

fp=[fpx,fpy];
fR=[fRx,fRy];

f=fp+fR;




end

