function [prf] = naca4(prf,nums)
%NACA4 erstellt Knoten des Profils
%   Speichert im profile-struct die Knoten Koordinaten prf.nodes.X und
%   prf.nodes.Y und die gesamte Knotenanzahl prf.M

c=prf.c;
m=nums(1)/100;  %Profilwoelbung
p=nums(2)/10;   %Woelbungsruecklage
t=(nums(3)*10+nums(4))/100; % Verhaeltnis Profildicke zu Sehnenlaenge

a=3;
x = 0.5*c*(tanh(a*(2*((0:(prf.N-1)))/(prf.N-1)-1))/tanh(a))+0.5*c;%linspace(0,c,prf.N);
yc = (x<=p).*(c*m/(p^2)*(2*p*x-x.^2)) + (x>p).*(c*m/((1-p)^2)*((1-2*p)+2*p*x-x.^2));
yt = 5*c*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);

x = x*c;
% mit Profilkrümmung
%theta = atan((x<=p).*(2*m/(p^2)*(p-x)) + (x>p).*(2*m/((1-p)^2)*(p-x))); 
% xu = x-yt.*sin(theta);             yu = yc+yt.*cos(theta);
% xl = x+yt.*sin(theta);             yl = yc-yt.*cos(theta);


% mit Profilkrümmung vernachlässigt
xu = x;             yu = yc+yt;
xl = x;             yl = yc-yt;


prf.nodes.X = [xu(end:-1:1) xl(2:end)];   prf.nodes.Y = [yu(end:-1:1) yl(2:end)];
prf.M = length(prf.nodes.X);

end