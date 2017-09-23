function p=create_panels(p)

p.panels.X=zeros(p.N,2);
p.panels.Y=zeros(p.N,2);
for i=1:p.M
    p.panels.X(i,:) = [p.nodes.X(i) p.nodes.X(mod(i,p.M)+1)]; 
    p.panels.Y(i,:) = [p.nodes.Y(i) p.nodes.Y(mod(i,p.M)+1)]; 
end