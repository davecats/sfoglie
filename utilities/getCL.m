function [CL , dCL_dalfa] = getCL(prf,U, AlfaDer)
%GETCL  calculates the lift coefficient for the profile with given boundary
%       edge tangential velocity using midpoint integration rule ( second order accurate)



h=( prf.panels.X(2,:)-prf.panels.X(1,:) )*cos(prf.alfa) + ( prf.panels.Y(2,:)-prf.panels.Y(1,:) )*sin(prf.alfa);
h=h';

SQU=(U(1:prf.N)/prf.Uinfty).^2;
SQU=[SQU;SQU(1)];

% mid point integral rule
CC1= 1- SQU(1:end-1) ;
CC2= 1- SQU(2:end  ) ;
CC= h.*(CC1+CC2)/2;

CL=sum(CC);


if nargin==3 && AlfaDer
    dxa=-(prf.panels.X(2,:)-prf.panels.X(1,:)  )*sin(prf.alfa) + (prf.panels.Y(2,:)-prf.panels.Y(1,:) )*cos(prf.alfa);
    CA= dxa'.*(CC1+CC2)/2;
    dCL_dalfa=sum(CA);
end

end

% without TE
% h=( prf.nodes.X(2:end)-prf.nodes.X(1:end-1) )*cos(prf.alfa) + ( prf.nodes.Y(2:end)-prf.nodes.Y(1:end-1) )*sin(prf.alfa);
% h=h';
% 
% SQU=(U(1:prf.N)/prf.Uinfty).^2;
% 
% % mid point integral rule
% CC1= 1- SQU(1:end-1) ;
% CC2= 1- SQU(2:end  ) ;
% CC= h.*(CC1+CC2)/2;
% 
% CL=sum(CC);
