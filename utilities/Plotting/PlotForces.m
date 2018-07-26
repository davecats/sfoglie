function [forces] = PlotForces(prf,flo,sol,mode, sol2)
%PLOTFORCES plots the line forces over the profile
%           mode=0: sum of lower and upper airfoil part
%           mode=1: upper airfoil part
%           mode=2: lower airfoil part
%           if sol2 comitted -> plot both solutions


if nargin<5
    Comp=false;
else
    Comp=true;
end


forces = LineForceDensity(prf,flo,sol);

if Comp
   forces2 = LineForceDensity(prf,flo,sol2); 
end

if mode==0 % sum of upper and lower part
    x=prf.nodes.X(1:prf.M);
    
    tmp=spline(prf.nodes.X(prf.M:prf.N),forces.fx_cum(prf.M:prf.N),x);
    xCum=forces.fx_cum(1:prf.M)+tmp';
    
    tmp=spline(prf.nodes.X(prf.M:prf.N),forces.fy_cum(prf.M:prf.N),x);
    yCum=forces.fy_cum(1:prf.M)+tmp';
    
    tmp=spline(prf.nodes.X(prf.M:prf.N),forces.fRx_cum(prf.M:prf.N),x);
    xRCum=forces.fRx_cum(1:prf.M)+tmp';
    
    tmp=spline(prf.nodes.X(prf.M:prf.N),forces.fRy_cum(prf.M:prf.N),x);
    yRCum=forces.fRy_cum(1:prf.M)+tmp';

    
    tmp=spline(prf.nodes.X(prf.M:prf.N),forces.fpx_cum(prf.M:prf.N),x);
    xpCum=forces.fpx_cum(1:prf.M)+tmp';
    
    tmp=spline(prf.nodes.X(prf.M:prf.N),forces.fpy_cum(prf.M:prf.N),x);
    ypCum=forces.fpy_cum(1:prf.M)+tmp';
    add='';
    
    if Comp
        tmp=spline(prf.nodes.X(prf.M:prf.N),forces2.fx_cum(prf.M:prf.N),x);
        xCum2=forces2.fx_cum(1:prf.M)+tmp';

        tmp=spline(prf.nodes.X(prf.M:prf.N),forces2.fy_cum(prf.M:prf.N),x);
        yCum2=forces2.fy_cum(1:prf.M)+tmp';

        tmp=spline(prf.nodes.X(prf.M:prf.N),forces2.fRx_cum(prf.M:prf.N),x);
        xRCum2=forces2.fRx_cum(1:prf.M)+tmp';

        tmp=spline(prf.nodes.X(prf.M:prf.N),forces2.fRy_cum(prf.M:prf.N),x);
        yRCum2=forces2.fRy_cum(1:prf.M)+tmp';


        tmp=spline(prf.nodes.X(prf.M:prf.N),forces2.fpx_cum(prf.M:prf.N),x);
        xpCum2=forces2.fpx_cum(1:prf.M)+tmp';

        tmp=spline(prf.nodes.X(prf.M:prf.N),forces2.fpy_cum(prf.M:prf.N),x);
        ypCum2=forces2.fpy_cum(1:prf.M)+tmp'; 
    end
    
   
elseif mode==1 % upper part
    x=prf.nodes.X(1:prf.M);
    
    xCum=forces.fx_cum(1:prf.M);
    yCum=forces.fy_cum(1:prf.M);
    xRCum=forces.fRx_cum(1:prf.M);
    yRCum=forces.fRy_cum(1:prf.M);
    xpCum=forces.fpx_cum(1:prf.M);
    ypCum=forces.fpy_cum(1:prf.M);
    add=' on upper part';
    
    if Comp
        xCum2=forces2.fx_cum(1:prf.M);
        yCum2=forces2.fy_cum(1:prf.M);
        xRCum2=forces2.fRx_cum(1:prf.M);
        yRCum2=forces2.fRy_cum(1:prf.M);
        xpCum2=forces2.fpx_cum(1:prf.M);
        ypCum2=forces2.fpy_cum(1:prf.M);  
    end 
else % lower part
    x=prf.nodes.X(prf.M:prf.N);
    
    xCum=forces.fx_cum(prf.M:prf.N);
    yCum=forces.fy_cum(prf.M:prf.N);
    xRCum=forces.fRx_cum(prf.M:prf.N);
    yRCum=forces.fRy_cum(prf.M:prf.N);
    xpCum=forces.fpx_cum(prf.M:prf.N);
    ypCum=forces.fpy_cum(prf.M:prf.N); 
    add=' on lower part';
    if Comp
        xCum2=forces2.fx_cum(prf.M:prf.N);
        yCum2=forces2.fy_cum(prf.M:prf.N);
        xRCum2=forces2.fRx_cum(prf.M:prf.N);
        yRCum2=forces2.fRy_cum(prf.M:prf.N);
        xpCum2=forces2.fpx_cum(prf.M:prf.N);
        ypCum2=forces2.fpy_cum(prf.M:prf.N);  
    end
    
    
    
end
str={'total forces','shear forces','pressure forces'};

figure
hold on
plot(x,xCum,'b')
plot(x,xRCum,'r')
plot(x,xpCum,'g')
if Comp
   plot(x,xCum2,'b --')
   plot(x,xRCum2,'r -- ')
   plot(x,xpCum2,'g --')
end
title(['cummulative streamwise forces',add])
xlabel('x')
ylabel('f_d_r_a_g')
legend(str,'location','best')

figure
hold on
plot(x,yCum,'b')
plot(x,yRCum,'r')
plot(x,ypCum,'g')
if Comp
   plot(x,yCum2,'b --')
   plot(x,yRCum2,'r --')
   plot(x,ypCum2,'g --') 
end
title(['cummulative lift forces',add])
xlabel('x')
ylabel('f_L')
legend(str,'location','best')


end

