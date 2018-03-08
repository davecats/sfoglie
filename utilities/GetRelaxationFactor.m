function Rel = GetRelaxationFactor( inp, rel0, minVal,maxVal )
%GETRELAXATIONFACTOR calculate factor for underrelaxation if needed


if nargin==1
    rel0=1;
elseif nargin==2
    maxVal=1.5;
    minVal=-0.5;
elseif nargin==3
    maxVal=1.5;
end
    

cur=rel0*inp;
Rel=ones(size(inp));
Rel(cur> maxVal)=  maxVal./inp(cur> maxVal);
Rel(cur< minVal)=  minVal./inp(cur< minVal);
%Rel=min(Rel);
Rel=min([Rel;rel0]);


end

