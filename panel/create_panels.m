function p=create_panels(p)
%CREATE_PANELS definiert die Randelemente
%   Speichert X und Y Koordinaten von Start- und Endpunkt des Randelements
%   profile.panel.X -> Dimension M x 2, 1. Spalte X-Koord Startpunkt, 2.
%   Spalte X-Koord Endpunkt (Y analog)


p.panels.X=zeros(p.N,2);
p.panels.Y=zeros(p.N,2);
p.panels.centre.X=zeros(p.N,1);
p.panels.length=zeros(p.N,1);
p.panels.centre.Y=zeros(p.N,1);
p.panels.tangent=zeros(p.N,2);
p.panels.normal=zeros(p.N,2);
for i=1:p.M
    p.panels.X(i,:) = [p.nodes.X(i) p.nodes.X(mod(i,p.M)+1)]; 
    p.panels.Y(i,:) = [p.nodes.Y(i) p.nodes.Y(mod(i,p.M)+1)]; 
    p.panels.centre.X(i,1) = mean([p.nodes.X(i) p.nodes.X(mod(i,p.M)+1)]);
    p.panels.centre.Y(i,1) = mean([p.nodes.Y(i) p.nodes.Y(mod(i,p.M)+1)]);
    t_vector  = [p.nodes.X(mod(i,p.M)+1)-p.nodes.X(i) p.nodes.Y(mod(i,p.M)+1)-p.nodes.Y(i) 0];
    n_vector = cross(t_vector, [0 0 1]);
    p.panels.length(i,1) = norm(t_vector);
    p.panels.tangent(i,:) = t_vector(1:2)./norm(t_vector(1:2));
    p.panels.normal(i,:) = n_vector(1:2)./norm(n_vector(1:2));
end