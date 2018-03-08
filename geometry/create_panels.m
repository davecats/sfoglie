function p=create_panels(p)
%CREATE_PANELS Create the panels for boundary element method. Adds the panel struct to profilestruct.
%              Saves each panels start and end point coordinates, the panel
%              length vector L, the arclength vector s, the panel angle
%              vector theta aswell as tangential vectors e and normalvectors n of each panel.
%              first line of p.panels.X gives X-coordinates of startpoints, second line the ones of the end point


%identify if TE is sharp or not
if p.nodes.X(1)==p.nodes.X(end) && p.nodes.Y(1)==p.nodes.Y(end)
    p.IsSharp=true;
    x1=p.nodes.X(1:end-1);
    x2=p.nodes.X(2:end);
    y1=p.nodes.Y(1:end-1);
    y2=p.nodes.Y(2:end);  
    p.N=p.N-1;
else
    p.IsSharp=false;
    x1=p.nodes.X(1:end);
    x2=[p.nodes.X(2:end), p.nodes.X(1)];
    y1=p.nodes.Y(1:end);
    y2=[p.nodes.Y(2:end), p.nodes.Y(1)];   
end


p.panels.X = [x1; x2]; 
p.panels.Y = [y1; y2]; 

e= [x2-x1; y2-y1];
% panel lengths
L= ( e(1,:).^2 + e(2,:).^2 ).^0.5;
% tangential vectors
e=[e(1,:)./L; e(2,:)./L];
% normal vectors
n=[-e(2,:);e(1,:)];
% arc length
s=[0,cumsum(L(1:end-1))];

% angle of panel tangent vector to x- axis

%p.panels.phi=phi2pi(e);
% theta(phi>pi/2)=phi(phi>pi/2)-3/2*pi;
% theta(phi<=pi/2)=phi(phi<=pi/2)+1/2*pi;
% theta=transpose(theta);
% theta=atan2(e(:,1),-e(:,2));
p.panels.theta=atan2(e(1,:),-e(2,:));

if p.IsSharp
    p.panels.theta(end)=pi;
    p.SdotT=0;
    p.ScrossT=1;
    p.panels.nn=(n(:,1:end-1)+n(:,2:end))/2;
    p.panels.nn=[n(:,1),p.panels.nn, n(:,end)];
else
    eTE=e(:,end); % trailing edge tangent vector
    e1=e(:,1); % tangent vector first panel
    eN=e(:,end-1); % tangent vector last panel
    sE=0.5*(eN-e1);
    p.SdotT = norm(eTE.*sE); 
    %s_cross_t = sqrt(1-s_t^2);
    p.ScrossT=abs(sE(1)*eTE(2)-sE(2)*eTE(1));
    
    p.nodes.n=(n(:,1:end-2)+n(:,2:end-1))/2;
    p.nodes.n=[n(:,1),p.nodes.n, n(:,end-1)];
end


p.panels.n=n;
p.panels.e=e;
p.s=s;
p.panels.L=L;

% TE Gap
p.gap=p.ScrossT*p.panels.L(end);


end

