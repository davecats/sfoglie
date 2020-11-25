function [prf] = geometry(prf)

%% load geometry
geo = load(prf.geometryFilename);
c=prf.c;        % chord

%% extract top and bottom surface
[~, i0] = min(geo(:,1));
top = geo(1:i0,:);
bottom = geo(i0:end,:);

%% get x distribution for node spread
if prf.pmode==1
    a=1.5;
    ap=a+1;

    i= 1:prf.M;
    ratio= (i-1)/(prf.M-1);
    x=1 - ap*ratio.*(1-ratio).^a-(1-ratio).^ap;
    x(end)=1;
    x=c*x;
elseif prf.pmode ==2
    a=3;
    x = 0.5*c*(tanh(a*(2*((0:(prf.M-1)))/(prf.M-1)-1))/tanh(a))+0.5*c;
end


%% interpolate nodes
yu = interp1(top(:,1),top(:,2),x,'pchip','extrap');
yl = interp1(bottom(:,1),bottom(:,2),x,'pchip','extrap');
prf.nodes.X = [x(end:-1:1) x(2:end)];   prf.nodes.Y = [yu(end:-1:1) yl(2:end)];
prf.N = length(prf.nodes.X);

%% check if profile has sharp TE
if norm(geo(1,:)-geo(end,:))<1e-6
    prf.sharpTE=true;
end
