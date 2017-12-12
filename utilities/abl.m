function [ erg ] = abl( u,x)
%ABL    calculates the derivate of u with central differences
%       x can be the vector of absolut coordinates of each node (if length(x)=length(u))
%       x can be the vector with all lengths of each subintervall (if length(x)=length(u)-1)

if length(x)==length(u) % x is the vector of coordinates
    %central difference
    erg=(u(3:1:end)-u(1:1:end-2)) ./ (x(3:1:end)-x(1:1:end-2)); 
    % left boundary (first node)
    erg=[(u(2)-u(1))/(x(2)-x(1)),erg]; % forward difference
    % right boundary (last node)
    erg=[erg,(u(end)-u(end-1))/(x(end)-x(end-1))]; % backward difference
    
elseif  length(x)==length(u)-1 % x is the vector of subintervall lengths
    %central difference
    erg=(u(3:1:end)-u(1:1:end-2)) ./ (x(1:1:end-1)+x(2:1:end)); 
 % left boundary (first node)
    erg=[(u(2)-u(1))/x(1),erg]; % forward difference
 % right boundary (last node)
    erg=[erg,(u(end)-u(end-1))/x(end)]; % backward difference

    % complete forward diff
%     erg=(u(2:1:end)-u(1:1:end-1)) ./ (x(1:1:end));     
%     % right boundary (last node)
%     erg=[erg,(u(end)-u(end-1))/x(end)]; % backward difference
end
end

