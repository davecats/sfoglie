function [prf] = naca4(prf,nums,WithoutCurvature)
%NACA4 calculates the nodes of NACA profiles. Creates profile struct with node X and Y coordinates
%      WithoutCurvature=true: the curvature of the profile is neglected

c=prf.c;
m=nums(1)/100;  %Profilwölbung
p=nums(2)/10;   %Wölbungsruecklage
t=(nums(3)*10+nums(4))/100; % Verhältnis Profildicke zu Sehnenlänge

a=3;
x = 0.5*c*(tanh(a*(2*((0:(prf.M-1)))/(prf.M-1)-1))/tanh(a))+0.5*c;%linspace(0,c,prf.M);
yc = (x<=p).*(c*m/(p^2)*(2*p*x-x.^2)) + (x>p).*(c*m/((1-p)^2)*((1-2*p)+2*p*x-x.^2));
yt = 5*c*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);

x = x*c;

if WithoutCurvature % Profilkrümmung vernachlässigt
    xu = x;             yu = yc+yt;
    xl = x;             yl = yc-yt;
else
    % mit Profilkrümmung
    theta = atan((x<=p).*(2*m/(p^2)*(p-x)) + (x>p).*(2*m/((1-p)^2)*(p-x))); 
    xu = x-yt.*sin(theta);             yu = yc+yt.*cos(theta);
    xl = x+yt.*sin(theta);             yl = yc-yt.*cos(theta);
end








prf.nodes.X = [xu(end:-1:1) xl(2:end)];   prf.nodes.Y = [yu(end:-1:1) yl(2:end)];
prf.N = length(prf.nodes.X);

end