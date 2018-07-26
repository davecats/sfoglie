function [ out] = LineForceDensity(prf,flo,sol)
%LINEFORCEDENSITY calculates the line force density over the panel in streamwise direction and perpendicular to it
%                 out.fL    : lift forces at each node
%                 out.fdrag : drag forces at each node
%                 out.fy_cum: cumulative lift forces over lower and upper part
%                 out.fx_cum: cumulative drag forces over lower and upper part
%                 out.fxR_cum, out.fyR_cum: shear contribution
%                 out.fxp_cum, out.fyp_cum: pressure contribution


% calculate stresses
[~,fR,fp]= StressVector(sol.Cp(1:prf.N),sol.tau(1:prf.N),prf.nodes.n',sol.Vb(1:prf.N) );

% Transform to coordinate system in stream direction
fRx=  fR(:,1)*cos(flo.alfa) + fR(:,2)*sin(flo.alfa);
fRy= -fR(:,1)*sin(flo.alfa) + fR(:,2)*cos(flo.alfa);

fpx=  fp(:,1)*cos(flo.alfa) + fp(:,2)*sin(flo.alfa);
fpy= -fp(:,1)*sin(flo.alfa) + fp(:,2)*cos(flo.alfa);

out.fL   = fRy+fpy;
out.fdrag= fRx+fpx;

% take half length of each panel adjanced to the node
L=0.5*(prf.panels.L(1:end-2)+prf.panels.L(2:end-1) );
% first and last node
L= [ prf.panels.L(1)/2; L' ; prf.panels.L(end-1)/2 ];

% multiply with panel lengths
out.fL   =out.fL   .*L;
out.fdrag=out.fdrag.*L;

% arclength vectors of lower and upper part starting at leading edge
sU= [0,cumsum(prf.panels.L(prf.M-1:-1:1))];
sL= [0,cumsum(prf.panels.L(prf.M:end-1))];

indU=prf.M:-1:1;
indL=prf.M:prf.N;


[~,t1]=NumInt(fRx(indU),sU,1,true);
[~,t2]=NumInt(fRx(indL),sL,1,true);
out.fRx_cum= [t1(indU);t2];

[~,t1]=NumInt(fRy(indU),sU,1,true);
[~,t2]=NumInt(fRy(indL),sL,1,true);
out.fRy_cum= [t1(indU);t2];

[~,t1]=NumInt(fpx(indU),sU,1,true);
[~,t2]=NumInt(fpx(indL),sL,1,true);
out.fpx_cum= [t1(indU);t2];

[~,t1]=NumInt(fpy(indU),sU,1,true);
[~,t2]=NumInt(fpy(indL),sL,1,true);
out.fpy_cum= [t1(indU);t2];

out.fx_cum=out.fRx_cum + out.fpx_cum;
out.fy_cum=out.fRy_cum + out.fpy_cum;

end

