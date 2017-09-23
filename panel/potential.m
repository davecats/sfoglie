function fld = potential(prf)

nodes = [mean(prf.panels.X,2) mean(prf.panels.Y,2)];  M = length(nodes(:,1));
A = zeros(M+1,M+1); A(1:M-1,M+1)=-1; t=zeros(M+1,1);

for i = 1:M
    x = nodes(i,1); y=nodes(i,2); 
    for j=1:M
        e=[prf.panels.X(j,2) prf.panels.Y(j,2)]-[prf.panels.X(j,1) prf.panels.Y(j,1)]; e=e/norm(e);
        r1=norm([prf.panels.X(j,1) prf.panels.Y(j,1)]-nodes(i));
        r2=norm([prf.panels.X(j,2) prf.panels.Y(j,2)]-nodes(i));
        x1=sum((nodes(i)-[prf.panels.X(j,1) prf.panels.Y(j,1)]).*e);
        x2=sum((nodes(i)-[prf.panels.X(j,2) prf.panels.Y(j,2)]).*e);
        y=cross([nodes(i,1) nodes(i,2) 0]-[prf.panels.X(j,2) prf.panels.Y(j,2) 0],[e 0]); y=y(3);
        t2=atan2(x-prf.panels.X(j,2),y-prf.panels.Y(j,2));
        t1=atan2(x-prf.panels.X(j,1),y-prf.panels.Y(j,1));
        psijp=x1*log(r1)-x2*log(r2)+x2-x1+y*(t1-t2);
        psijm=( (x1+x2)*psijp + r2^2*log(r2) -r1^2*log(r1)+0.5*(x1^2-x2^2) )/(x1-x2);
        %psis
        A(i,[j j+1])=A(i,[j j+1])+psijp/(4*pi); 
        A(i,[j j+1])=A(i,[j j+1])+[-1 1]*psijm/(4*pi);
    end
    t(i)=-1*y;
end
e=[prf.panels.X(end,2) prf.panels.Y(end,2)]-[prf.panels.X(end,1) prf.panels.Y(end,1)]; e=e/norm(e);
t1=atan2(prf.panels.Y(1,2)-prf.panels.Y(1,1),prf.panels.X(1,2)-prf.panels.X(1,1));
t2=atan2(prf.panels.Y(end,1)-prf.panels.Y(end,2),prf.panels.X(end,1)-prf.panels.X(end,2)); t2=t2+(t2<0)*2*pi;
s=[cos(t2) sin(t2)];
A(M,[1 M-1])=0.5*([1 -1])*norm(e.*s);
A(end,[1, end-1])=[1,1];


fld.gamma = A\t;
fld.nodes = nodes;

end