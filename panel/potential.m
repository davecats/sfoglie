function fld = potential(prf)

nodes = [mean(prf.panels.X,2) mean(prf.panels.Y,2)];  
M = length(nodes(:,1));    % number of panels (including trailing edge)
N = M-1;                   % number of nodes
A = zeros(N+1,N+1); A(1:N,N+1)=-1; 
t=zeros(N+1,1);

% Unit vector of last panel
e=[prf.panels.X(end,2) prf.panels.Y(end,2)]-[prf.panels.X(end,1) prf.panels.Y(end,1)]; e=e/norm(e);
t1=atan2(prf.panels.Y(1,2)-prf.panels.Y(1,1),prf.panels.X(1,2)-prf.panels.X(1,1));                 t1=t1+(t1<0)*2*pi;
t2=atan2(prf.panels.Y(end-1,1)-prf.panels.Y(end-1,2),prf.panels.X(end-1,1)-prf.panels.X(end-1,2)); t2=t2+(t2<0)*2*pi;
s=-[cos(0.5*(t2+t1)) sin(0.5*(t2+t1))];
s_t = norm(e.*s);
% Influence matrix for i=M-1 control nodes and j=M panels
for i = 1:N
    x = nodes(i,1); y=nodes(i,2); 
    %x=prf.panels.X(i,1); y=prf.panels.Y(i,1);
    for j=1:M
        e=[prf.panels.X(j,2) prf.panels.Y(j,2)]-[prf.panels.X(j,1) prf.panels.Y(j,1)]; e=e/norm(e);
        n=[-e(2) e(1)];
        r1=norm([prf.panels.X(j,1) prf.panels.Y(j,1)]-nodes(i,:));
        r2=norm([prf.panels.X(j,2) prf.panels.Y(j,2)]-nodes(i,:));
        x1=-sum((nodes(i,:)-[prf.panels.X(j,1) prf.panels.Y(j,1)]).*e);
        x2=-sum((nodes(i,:)-[prf.panels.X(j,2) prf.panels.Y(j,2)]).*e);
        Y=-sum((nodes(i,:)-[prf.panels.X(j,1) prf.panels.Y(j,1)]).*n);
        
        t1=atan2(x-prf.panels.X(j,1),y-prf.panels.Y(j,1));  t1=t1+(t1<0)*2*pi; 
        t2=atan2(x-prf.panels.X(j,2),y-prf.panels.Y(j,2));  t2=t2+(t2<0)*2*pi; if t1*t2<0; if t1<0; t1=t1+2*pi; else t2=t2+2*pi; end; end;
        t1_t2=t1-t2; %while t1_t2<-2*pi; t1_t2=t1_t2-2*pi; end; %fld.t1pt2(i,j)=t1+t2; fld.t1(i,j)=t1; fld.t2(i,j)=t2; fld.t1mt2(i,j)=t1-t2;
        
        %t1=atan2(y-prf.panels.Y(j,1),x-prf.panels.X(j,1));  t1=t1+(t1<0)*2*pi; t1=pi/2-t1; t1=t1+(t1<0)*2*pi; t1=t1-(t1>pi)*2*pi;
        %t2=atan2(y-prf.panels.Y(j,2),x-prf.panels.X(j,2));  t2=t2+(t2<0)*2*pi; t2=pi/2-t2; t2=t2+(t2<0)*2*pi; t2=t2-(t2>pi)*2*pi; 
        %t1_t2=t1-t2; %while t1_t2>2*pi; t1_t2=t1_t2-2*pi; end; 
        psijp=x1*log(r1)-x2*log(r2)+x2-x1+Y*(t1_t2);
        psijm=( (x1+x2)*psijp + (r2^2)*log(r2) -(r1^2)*log(r1)+0.5*(x1^2-x2^2) )/(x1-x2);
        if j<N
          %psi+
          A(i,[j j+1]) = A(i,[j j+1]) +psijp*[ 1 1]/(4*pi); 
          %psi-
          A(i,[j j+1]) = A(i,[j j+1]) +psijm*[-1 1]/(4*pi);
        else
          A(i,[1 N])   = A(i,[1 N])   +psijp*s_t*[1 -1]/(4*pi); 
        end
    end
    t(i)=-1*y;
end
A(N+1,[1, N])=[1,1];

fld.A=A;
fld.t=t;
fld.gamma = A\t;
fld.nodes = nodes;

end