function fld = potential(prf)

nodes = [mean(prf.panels.X,2) mean(prf.panels.Y,2)];  
M = length(nodes(:,1));    % number of panels (including trailing edge)
N = M-1;                   % number of panels (excluding trailing edge)
A = zeros(N+2,N+2); A(1:N+1,N+2)=-1; 
t=zeros(N+2,1);
u = cos(prf.alfa);
v = sin(prf.alfa);

% Unit vector of last panel
e=[prf.panels.X(end,2) prf.panels.Y(end,2)]-[prf.panels.X(end,1) prf.panels.Y(end,1)]; e=e/norm(e);
if (atan2(e(2),e(1))>pi/2 || atan2(e(2),e(1))<-pi/2); e=-e; end
e1=[prf.panels.X(1,2) prf.panels.Y(1,2)]-[prf.panels.X(1,1) prf.panels.Y(1,1)];  e1=e1/norm(e1);
eN=-[prf.panels.X(N,2) prf.panels.Y(N,2)]+[prf.panels.X(N,1) prf.panels.Y(N,1)]; eN=eN/norm(eN);
s=0.5*(e1+eN); s=s/norm(s); s_t = norm(e.*s); s_cross_v = sqrt(1-s_t^2);

% Influence matrix for i=M-1 control nodes and j=M panels
for i = 1:N+1
    %x = nodes(i,1); y=nodes(i,2); 
    x = prf.panels.X(i,1); y=prf.panels.Y(i,1);
    for j=1:N+1
         [psi1,psi2]=linevortex([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                                 prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                                [x,y]);    
        if j<=N
          A(i,[j j+1]) = A(i,[j j+1]) + [psi1 psi2]; 
        else
          [psi]=linesource([prf.panels.X(j,1) prf.panels.Y(j,1);  ...
                            prf.panels.X(j,2) prf.panels.Y(j,2)], ... 
                           [x,y]);
          A(i,[1 N+1]) = A(i,[1 N+1]) + 0.5*s_t*[1 -1]*(psi1+psi2) ...
                                      + 0.5*s_cross_v*[1 -1]*psi;         
        end
    end
    t(i)=-u*y+v*x;
end
A(N+2,[1, N+1])=[1,1];

fld.A=A;
fld.t=t;
fld.gamma = A\t;
fld.nodes = nodes;

end