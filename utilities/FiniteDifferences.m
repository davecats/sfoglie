function [ dV ] = FiniteDifferences( V,x, verf )
%FINITEDIFFERENCES calculates the finite difference approximation of a vector. 
%                       verf='central'  -> uses central differences 
%                       verf='backward' -> uses backward differences 
%                       verf='forward'  -> uses central differences 
%                   first vector element always calculated with forward differences
%                   last vector element always calculated with backward differencest


if nargin==2
    verf='central';
end

if length(V(:,1))==1 && length(x(:,1))==1
    V=transpose(V);x=transpose(x);
elseif length(V(:,1))==1
    V=transpose(V);
elseif length(x(:,1))==1
    x=transpose(x);
elseif length(V)~=length(x)
    disp('FiniteDifferences: dimensions of vectors do not fit')
    return
end

if strcmp(verf,'forward')
    dV= (V(2:end)-V(1:end-1)) ./ (x(2:end)-x(1:end-1) );
    dV=[dV ; dV(end)]; % last node backward difference
elseif strcmp(verf,'backward')
    dV= (V(2:end)-V(1:end-1)) ./ (x(2:end)-x(1:end-1) );   
    dV=[dV ; dV(end)]; % first node forward difference
else % central differences
    dV= (V(3:end)-V(1:end-2)) ./ (x(3:end)-x(1:end-2) );
    dV1= ( V(2)-V(1) )/(x(2)-x(1)); % first node forward difference
    dVN= ( V(end)-V(end-1) )/(x(end)-x(end-1)); % last node backward difference
    dV= [dV1; dV; dVN];
end


end

