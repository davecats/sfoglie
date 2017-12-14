function [fld] = potential(prf)
%POTENTIAL  calculates the coefficient matrix A of the equation system 
%           A_ij z_j = -U y_i + V x_i with solution vector z=[gamma_1,..gamma_M, psi_0]
%           and solves it


%nodes = [mean(prf.panels.X,2) mean(prf.panels.Y,2)];  %linewise mean value -> panel mid point
%fld.nodes = nodes;
M = length(prf.panels.X);    % number of panels (including trailing edge)
N = M-1;                   % number of panels (excluding trailing edge)

B = zeros(N+1,N+1); % Matrix for q_j -> psi_q= B q 

A = zeros(N+2,N+2);  % Matrix for gamma_j -> psi_gam= A g 
A(1:N+1,N+2)=-1; % element N+2 of the solution vector is the unknown psi_0 
                 %  ->appears in every equation besides the Kutta condition
                             
t=zeros(N+2,1); % right hand side vector
u=evalin('base','ui');
v= evalin('base','vi');



eTE=prf.panels.e(end,:); % trailing edge tangent vector
e1=prf.panels.e(1,:); % tangent vector first panel
eN=prf.panels.e(end-1,:); % tangent vector last panel
s=0.5*(eN-e1);
s_t = norm(eTE.*s); s_cross_v = sqrt(1-s_t^2);

fld.SdotT=s_t;
fld.ScrossT=s_cross_v;

TE= zeros(N+2,1);%Trailing edge part of equation System



% Influence matrix for i=M control nodes and j=M panels (including trailing edge panel)
for i = 1:N+1 % going through all nodes
    %x = nodes(i,1); y=nodes(i,2); 
    x = prf.panels.X(i,1); y=prf.panels.Y(i,1); % Kollocation method: loading point in i-th node -> i-th equation
    for j=1:N+1 % sum over all Panels ->coefficient of gamma_j in i-th eqation
         % Coefficients of gamma_j and gamma_j+1 from j-th panel -> linear ansatz 
         [psi1,psi2]=linevortex([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                                 prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                                [x,y]); 
         % Coefficients of q_j from j-th panel -> constant ansatz                  
         [psi]=linesource([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                            prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                           [x,y]);        
        if j<=N 
          A(i,[j j+1]) = A(i,[j j+1]) + [psi1 psi2]; 
        else %trailing edge panel N+1
            TE(i)= s_t*(psi1+psi2)-s_cross_v*psi;
          %A(i,1) = A(i,1) + s_t*(psi1+psi2) ...  % psi^+= psi_1 + psi_2 
          %                - s_cross_v*psi;       
        end
        B(i,j)=psi;
    end
    t(i)=-u*y+v*x; % free stream -> right hand side
end
A(N+2,[1, N+1])=[1,1];%Kutta condition
A(:,1)=A(:,1)+TE(:);

fld.A=A(1:N+1,1:N+1); % without the psi0 part
fld.B=B;
fld.t=t(1:N+1);
Lsg=A\t;    %solve the equation system
fld.gamma = Lsg(1:N+1); % gammas -> solution vector without psi0 -> last element
fld.psi0=Lsg(N+2);


end