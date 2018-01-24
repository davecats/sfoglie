function [X1,X2,Y] =GetLocalRelCoord(xi,eta,prf,withoutTE)
%EVALUATEBOUNDARYINT calculates the relativ coordinates in local panel coordinate system for airfoil panels
%                    xi:  column vector with Loadingpoint x-coordinate
%                    eta: column vector with Loadingpoint y-coordinate
%                    withoutTE: if comitted the trailing edge part is left

n=prf.panels.n; % normal vector of each panel
e=prf.panels.e; % tangent vector of each panel
N=prf.N; % number of nodes

M=length(xi);
% Panel node 1 and 2 coordinates
xp1=prf.panels.X(1,:);% 1- start point of panel
yp1=prf.panels.Y(1,:);
xp2=prf.panels.X(2,:);% 2- end point of panel
yp2=prf.panels.Y(2,:);

% relativ vectors in global KOS
dx1= xi*ones(1,N)  - ones(M,1)*xp1;
dx2= xi*ones(1,N)  - ones(M,1)*xp2;
dy1= eta*ones(1,N) - ones(M,1)*yp1;
dy2= eta*ones(1,N) - ones(M,1)*yp2;

% relativ vectors in local KOS
X1=(ones(M,1)*e(1,:)).*dx1+(ones(M,1)*e(2,:)).*dy1;
X2=(ones(M,1)*e(1,:)).*dx2+(ones(M,1)*e(2,:)).*dy2;
Y =(ones(M,1)*n(1,:)).*dx1+(ones(M,1)*n(2,:)).*dy1;

% no TE part
if nargin==4
    X1=X1(:,1:end-1);
    X2=X2(:,1:end-1);
    Y=Y(:,1:end-1);
end

end

