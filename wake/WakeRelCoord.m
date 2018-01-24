function [X1,X2,Y] = WakeRelCoord( xi,eta,wake )
%WAKERELCOORD   calculates the relativ coordinates in local panel coordinate system for wake panels
%                    xi:  column vector with Loadingpoint x-coordinate
%                    eta: column vector with Loadingpoint y-coordinate
%                    withoutTE: if comitted the trailing edge part is left


n=wake.n; % normal vector of each panel
e=wake.e; % tangent vector of each panel

M=length(xi);
% Panel node 1 and 2 coordinates
xp1=transpose(wake.x(1:end-1));% 1- start point of panel
yp1=transpose(wake.y(1:end-1));
xp2=[xp1(2:end),wake.x(end)];  % 2- end point of panel
yp2=[yp1(2:end),wake.y(end)];

N=length(xp1);% number of wake panels

% relativ vectors in global KOS
dx1= xi*ones(1,N)  - ones(M,1)*xp1;
dx2= xi*ones(1,N)  - ones(M,1)*xp2;
dy1= eta*ones(1,N) - ones(M,1)*yp1;
dy2= eta*ones(1,N) - ones(M,1)*yp2;

% relativ vectors in local KOS
X1=(ones(M,1)*e(1,:)).*dx1+(ones(M,1)*e(2,:)).*dy1;
X2=(ones(M,1)*e(1,:)).*dx2+(ones(M,1)*e(2,:)).*dy2;
Y =(ones(M,1)*n(1,:)).*dx1+(ones(M,1)*n(2,:)).*dy1;

end

