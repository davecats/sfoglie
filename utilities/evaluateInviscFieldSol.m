function [ psi ] = evaluateInviscFieldSol(xi,fld,prf)
%EVALUATEINVISCFIELDSOL     evaluates the streamfunktion for an arbitrary
%                           Field point

M = length(fld.gamma);    % number of panels (including trailing edge)
N = M-1;  

A=zeros(1,N+1);
for j=1:N 
     % Coefficients of gamma_j and gamma_j+1 from j-th panel -> linear ansatz 
     [psi1,psi2]=linevortex([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                             prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                            [xi(1), xi(2)]); 
      A([j j+1]) = A([j j+1]) + [psi1 psi2];       
end

%insert TE contribution
[psiT1,psiT2]=linevortex([prf.panels.X(M,1) prf.panels.Y(M,1);  ...
                             prf.panels.X(M,2) prf.panels.Y(M,2)], ... 
                            [xi(1), xi(2)]); 
psiTE= psiT1+psiT2;
A(1)=A(1)-fld.SdotT*psiTE;
psiS=linesource([prf.panels.X(M,1) prf.panels.Y(M,1);  ...
                             prf.panels.X(M,2) prf.panels.Y(M,2)], ... 
                            [xi(1), xi(2)]); 
A(1)=A(1)+fld.ScrossT*psiS ;




psi=A*fld.gamma;
u=evalin('base','ui');
v= evalin('base','vi');
% add freee stream
psi=psi+u*xi(2)-v*xi(1);


end

