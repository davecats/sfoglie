function fld = potential(prf)
%POTENTIAL  calculates the coefficient matrix A of the equation system 
%           A_ij z_j = -U y_i + V x_i with solution vector z=[gamma_1,..gamma_M, psi_0]
%           and solves it


nodes = [mean(prf.panels.X,2) mean(prf.panels.Y,2)];  %linewise mean value -> panel mid point
M = length(nodes(:,1));    % number of panels (including trailing edge)
N = M-1;                   % number of panels (excluding trailing edge)
A = zeros(N+2,N+2); A(1:N+1,N+2)=-1; % element N+2 of the solution vector is the unknown psi_0 
                                     %->appears in evere equation besides the Kutta condition
t=zeros(N+2,1); % right hand side vector
u = cos(prf.alfa);
v = sin(prf.alfa);

% Unit vector of last panel
e=[prf.panels.X(end,2) prf.panels.Y(end,2)]-[prf.panels.X(end,1) prf.panels.Y(end,1)]; e=e/norm(e);
if (atan2(e(2),e(1))>pi/2 || atan2(e(2),e(1))<-pi/2); e=-e; end 
e1=[prf.panels.X(1,2) prf.panels.Y(1,2)]-[prf.panels.X(1,1) prf.panels.Y(1,1)];  e1=e1/norm(e1);
eN=-[prf.panels.X(N,2) prf.panels.Y(N,2)]+[prf.panels.X(N,1) prf.panels.Y(N,1)]; eN=eN/norm(eN);
s=0.5*(e1+eN); s=s/norm(s); s_t = norm(e.*s); s_cross_v = sqrt(1-s_t^2);

% Influence matrix for i=M control nodes and j=M panels (including trailing edge panel)
for i = 1:N+1 % going through all nodes
    %x = nodes(i,1); y=nodes(i,2); 
    x = prf.panels.X(i,1); y=prf.panels.Y(i,1); % Kollocation method: loading point in i-th node -> i-th equation
    for j=1:N+1 % sum over all Panels ->coefficient of gamma_j in i-th eqation
         [psi1,psi2]=linevortex([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                                 prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                                [x,y]);    
        if j<=N % trailing edge panel N+1
          A(i,[j j+1]) = A(i,[j j+1]) + [psi1 psi2]; 
        else
          [psi]=linesource([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                            prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                           [x,y]);
          A(i,1) = A(i,1) + s_t*(psi1+psi2) ...  % psi^+= psi_1 + psi_2 
                          - s_cross_v*psi;       
        end
    end
    t(i)=-u*y+v*x; % free stream -> right hand side
end
A(N+2,[1, N+1])=[1,1];%Kutta condition

fld.A=A;
fld.t=t;
fld.gamma = A\t; %solve the equation system
fld.nodes = nodes;

end