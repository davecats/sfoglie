function p=create_panels(p)
%CREATE_PANELS  saves length L, normal vector n, tangential vector e and start and
%               end point of each panel. profile.panel.X/Y is M x 2 matrix.
%               first column -> start node coordinates, second col.-> end node.
%               saves also the arclength of each node in vector s

x1=p.nodes.X(1:end);
x2=[p.nodes.X(2:end), p.nodes.X(1)];
y1=p.nodes.Y(1:end);
y2=[p.nodes.Y(2:end), p.nodes.Y(1)];

p.panels.X = [transpose(x1), transpose(x2)]; 
p.panels.Y = [transpose(y1), transpose(y2)]; 

e= [transpose(x2-x1), transpose(y2-y1)];
% panel lengths
L= ( e(:,1).^2 + e(:,2).^2 ).^0.5;
% tangential vectors
e=[e(:,1)./L, e(:,2)./L];
% normal vectors
n=[e(:,2),-e(:,1)];
% arc length
s=[0;cumsum(L(1:end-1))];


eTE=e(end,:); % trailing edge tangent vector
e1=e(1,:); % tangent vector first panel
eN=e(end-1,:); % tangent vector last panel
sE=0.5*(eN-e1);
s_t = norm(eTE.*sE); s_cross_t = sqrt(1-s_t^2);

p.SdotT=s_t;
p.ScrossT=s_cross_t;

p.panels.n=n;
p.panels.e=e;
p.panels.s=s;
p.panels.L=L;


%  old with loop

% for i=1:p.M
%     x1=p.nodes.X(i);
%     x2=p.nodes.X(mod(i,p.M)+1);
%     y1=p.nodes.Y(i);
%     y2=p.nodes.Y(mod(i,p.M)+1);
%     p.panels.X(i,:) = [x1, x2 ]; 
%     p.panels.Y(i,:) = [y1, y2]; 
%     
%     % panels tangential and normal vectors
%     
%     l=norm(e);
%     L(i)=l;
%     e=e/l;
%     p.panels.e(i,:)=e;
%     n=[e(2) -e(1)]; 
%     p.panels.n(i,:)=n;
%     s(i+1)=s(i)+l;
% end
%s=s(1:end-1);
%p.panels.s=s;
%p.panels.L=L;

