function [ ] = PlotStressVectors( prf,sol,flo,each_nth,total, pressure, shear, force )
%PLOTSTRESSVECTORS plots the stress vectors at the profile
%                  each_nth: only plot each n-th node (default 1 -> every node)
%                  total   =true: plot the total stresses consisting of pressure and shear stresses 
%                  pressure=true: plot the pressure stresses
%                  shear   =true: plot the shear stresses
%                  force = true : weight the stresses with panel length (default false) 

if nargin<8
   force=false ;
end
if nargin<7
    shear=false;
end
if nargin<6
    pressure=false;
end
if nargin<5
    total=true;
end
if nargin<4
   each_nth=1; 
end

% calculate stresses
[f,fR,fp] = StressVector(sol.Cp(1:prf.N),2*sol.tau(1:prf.N),prf.nodes.n',sol.Vb/flo.Uinfty );

if force
    % take half length of each panel adjanced to the node
    L=0.5*(prf.panels.L(1:end-2)+prf.panels.L(2:end-1) );
    % first and last node
    L= [ prf.panels.L(1)/2; L' ; prf.panels.L(end-1)/2 ];
    L=[L,L];
    f=f.*L;
    fR=fR.*L;
    fp=fp.*L;
end


% scale for plot
if shear
    scale= 0.02*prf.c/max(max(fR));
end
if pressure
    scale= 0.2*prf.c/max(max(fp));
end
if total
    scale= 0.2*prf.c/max(max(f));
end

startX=prf.nodes.X';
startY=prf.nodes.Y';

legendstr={};
figure
hold on
plot([prf.panels.X],[prf.panels.Y],'k','Linewidth',1);

if total
    endX_tot= startX + scale*f(:,1);
    endY_tot= startY + scale*f(:,2);
    for i=1:each_nth:prf.N
        line([startX(i) endX_tot(i)]  , [startY(i) endY_tot(i)],'color','b');
    end
    %legendstr{end+1}='stress vector';
end


if shear
    endX_R= startX + scale*fR(:,1);
    endY_R= startY + scale*fR(:,2);
    for i=1:each_nth:prf.N
        line([startX(i) endX_R(i)]  , [startY(i) endY_R(i)],'color','r');
    end
    %legendstr{end+1}='shear stresses';
end


if pressure
    endX_p= startX + scale*fp(:,1);
    endY_p= startY + scale*fp(:,2);
    for i=1:each_nth:prf.N
        line([startX(i) endX_p(i)]  , [startY(i) endY_p(i)],'color','r');
    end
    %legendstr{end+1}='pressure stresses';
end
title(['Profile stress vectors'])
axis equal; xlabel('x'); ylabel('y') 
%legend(legendstr)

end

