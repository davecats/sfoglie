function [prf] = naca4(prf,nums,WithoutCurvature, sharp, mode)
%NACA4 calculates the nodes of NACA profiles. Creates profile struct with node X and Y coordinates
%      WithoutCurvature=true: the curvature of the profile is neglected
%      sharp = true         : adjust factor for NACA formular to get a sharp Trailing edge
%      mode=1 (default) / 2 : sets different node distributions over the chord length
%                           -> mode 1 more nodes in the middle of the profile
%                           -> mode 2 more nodes at Leading and Trailing edge, less in middle


if nargin<5
    mode=1;
end
if nargin<4
    sharp=false;
end
if nargin<3
    WithoutCurvature=false;
end

c=prf.c;        % chord
m=nums(1)/100;  % chamber
p=nums(2)/10;   % chamber position
t=(nums(3)*10+nums(4))/100; % ratio between maximum thickness and chord


% get x distribution for node spread
if mode==1
    a=1.5;
    ap=a+1;

    i= 1:prf.M;
    ratio= (i-1)/(prf.M-1);
    x=1 - ap*ratio.*(1-ratio).^a-(1-ratio).^ap;
    x(end)=1;
    x=c*x;
elseif mode ==2
    a=3;
    x = 0.5*c*(tanh(a*(2*((0:(prf.M-1)))/(prf.M-1)-1))/tanh(a))+0.5*c;
end
    
% adjusts coefficient of last Term in case of sharp Trailing edge
if sharp; coeff=0.1036; else coeff=0.1015; end
yt = 5*c*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - coeff*x.^4);

% additional chamber y-term
yc = (x<=p).*(c*m/(p^2)*(2*p*x-x.^2)) + (x>p).*(c*m/((1-p)^2)*((1-2*p)+2*p*x-x.^2));

x = x*c;
if WithoutCurvature % curvature neglected
    xu = x;             yu = yc+yt;
    xl = x;             yl = yc-yt;
else
    theta = atan((x<=p).*(2*m/(p^2)*(p-x)) + (x>p).*(2*m/((1-p)^2)*(p-x))); 
    xu = x-yt.*sin(theta);             yu = yc+yt.*cos(theta);
    xl = x+yt.*sin(theta);             yl = yc-yt.*cos(theta);
end



prf.nodes.X = [xu(end:-1:1) xl(2:end)];   prf.nodes.Y = [yu(end:-1:1) yl(2:end)];
prf.N = length(prf.nodes.X);


% Plot of node distributions for both modes
%-------------------------------------------------------------------------------

% m=80;
%     a=3;
%     x1 = 0.5*(tanh(a*(2*((0:(m-1)))/(m-1)-1))/tanh(a))+0.5;
% 
%     a=1.5;
%     ap=a+1;
% 
%     i= 1:m;
%     ratio= (i-1)/(m-1);
%     x2=1 - ap*ratio.*(1-ratio).^a-(1-ratio).^ap;
%     x2(end)=1;
% 
% 
% dx1=x1(2:end)-x1(1:end-1);
% dx2=x2(2:end)-x2(1:end-1);
% 
% x1m=0.5*(x1(2:end)+x1(1:end-1));
% x2m=0.5*(x2(2:end)+x2(1:end-1));
% 
% 
% figure
% hold on
% plot(x1m,dx1, 'b .')
% plot(x2m,dx2, 'g x')
% title('Boundary element length over midpoint x-koord')
% xlabel('x/c')
% ylabel('L_i')
% legend('mode 1','mode 2')

end