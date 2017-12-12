function [Cgam, Cq] = WakeNodeEqs(wake,prf)
% WAKENODEEQS   calculates wake velocity and Coefficients for the wake node
%               equations

N=evalin('base','N'); % number of nodes on airfoil
NW=evalin('base','NW');
%B=zeros(NW,N+NW);
Cgam=zeros(NW,N);
Cq=zeros(NW,N+NW);

for i=1:NW
    for j=1:N-1
        [c1, c2]=WakeCirc( [prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                           prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                          [wake.x(i),wake.y(i)],prf.panels.n(j,:));
        Cgam(i,j:j+1)=Cgam(i,j:j+1)+[c1, c2];
    end  
end


for i=2:NW
    [~,du_dq] = WakeSourceCoeffs( wake ,[wake.x(i),wake.y(i)] , wake.n(i-1,:));
    Cq(i,N+1:end)=du_dq;   
end




% for i=1:NW
%    for k=1:N-1 % constant source distribution for airfoil souerces
%        xi=[xw(i,1),xw(i,2)]; % loading point is i-th wake node
%        psi=linesource([prf.panels.X(k,1) prf.panels.Y(k,1);  ...
%                             prf.panels.X(k,2) prf.panels.Y(k,2)], ... 
%                            xi); 
%        B(i,k)=psi;
%    end
%     % linear source distribution for wake souerces
%    psi=GetWakeCoeffs( xw , xi);
%    B(i,N+1:end)=psi;
% end
% %clear xi psi;

end

