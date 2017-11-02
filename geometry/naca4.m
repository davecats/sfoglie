function [prf] = naca4(prf,nums)

c=prf.c;
m=nums(1)/100;
p=nums(2)/10;
t=(nums(3)*10+nums(4))/100;

a=3;
x = 0.5*c*(tanh(a*(2*((0:(prf.N-1)))/(prf.N-1)-1))/tanh(a))+0.5*c;%linspace(0,1,prf.N);
yc = (x<=p).*(c*m/(p^2)*(2*p*x-x.^2)) + (x>p).*(c*m/((1-p)^2)*((1-2*p)+2*p*x-x.^2));
yt = 5*c*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);
theta = atan((x<=p).*(2*m/(p^2)*(p-x)) + (x>p).*(2*m/((1-p)^2)*(p-x))); 

x = x*c;
xu = x-yt.*sin(theta);             yu = yc+yt.*cos(theta);
xl = x+yt.*sin(theta);             yl = yc-yt.*cos(theta);

prf.nodes.X = [xu(end:-1:1) xl(2:end)];   prf.nodes.Y = [yu(end:-1:1) yl(2:end)];
prf.M = length(prf.nodes.X);

end