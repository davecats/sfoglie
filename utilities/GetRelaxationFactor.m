function Rel = GetRelaxationFactor( inp, rel0, minVal,maxVal )
%GETRELAXATIONFACTOR calculate factor for underrelaxation if needed
%                       -> underrelaxation if relative changes are bigger than defined limits
%                    inp: relativ changes (z_k+1 - z_k) / z_k for k-th iteration
%                    minVal: lower limit (default -0.5)
%                    maxVal: upper limit (default  1.5)
%                    rel0:   already applied relaxation factor (default 1)


if nargin==1
    rel0=1;
    maxVal=1.5;
    minVal=-0.5;
elseif nargin==2
    maxVal=1.5;
    minVal=-0.5;
elseif nargin==3
    maxVal=1.5;
end
    
if rel0<0
    rel0=0.1;
end

cur=rel0*inp;
Rel=ones(size(inp));
Rel(cur> maxVal)=  maxVal./inp(cur> maxVal);
Rel(cur< minVal)=  minVal./inp(cur< minVal);
Rel=min(Rel);
Rel=min([Rel;rel0]);


end

