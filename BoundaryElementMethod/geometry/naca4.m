function [prf] = naca4(prf)
%NACA4 calculates the node coordinates of NACA profiles. Creates profile struct with node X and Y coordinates
%      WithoutCurvature=true: the curvature of the profile is neglected
%      sharp = true         : adjust factor for NACA formular to get a sharp Trailing edge
%      mode=1 (default) / 2 : sets different node distributions over the chord length
%                           -> mode 1 more nodes in the middle of the profile
%                           -> mode 2 more nodes at Leading and Trailing edge, less in middle

c=prf.c;        % chord
m=prf.naca(1)/100;  % chamber
p=prf.naca(2)/10;   % chamber position
t=(prf.naca(3)*10+prf.naca(4))/100; % ratio between maximum thickness and chord


% get x distribution for node spread
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
    
% adjusts coefficient of last Term in case of sharp Trailing edge
if prf.sharpTE; coeff=0.1036; else coeff=0.1015; end
% profile thicknes depending on x
yt = 5*c*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - coeff*x.^4);

% Skeleton line of profile
yc = (x<=p).*(c*m/(p^2)*(2*p*x-x.^2)) + (x>p).*(c*m/((1-p)^2)*((1-2*p)+2*p*x-x.^2));

x = x*c;

% add thicknes to skeleton line to get the profile

if prf.noSkew % curvature neglected
    xu = x;             yu = yc+yt;
    xl = x;             yl = yc-yt;
else % 
    theta = atan((x<=p).*(2*m/(p^2)*(p-x)) + (x>p).*(2*m/((1-p)^2)*(p-x))); 
    xu = x-yt.*sin(theta);             yu = yc+yt.*cos(theta);
    xl = x+yt.*sin(theta);             yl = yc-yt.*cos(theta);
end

% brings nodes in counter clock wise order starting at the TE on the upper
% surface

prf.nodes.X = [xu(end:-1:1) xl(2:end)];   prf.nodes.Y = [yu(end:-1:1) yl(2:end)];
prf.N = length(prf.nodes.X);


end