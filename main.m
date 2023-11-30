%% Initialisation
close all
clear all

addpath('./panel/')
addpath('./geometry/')
addpath('./utilities/')
set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');


%% Parameters
profile.c = 1;
profile.N = 500; %(Oben und unten N-1 Panels + Hinterseite) -> M=1+2(N-1)=2N-1 Randelemente
profile.alfa = 2*pi/180;
profile.noSkew = true;
NACA = [4 4 1 2];

%% Create profile and panels
profile = naca4(profile,NACA);
profile = create_panels(profile);

%% Solve potential problem
field=potential(profile); Nle = round((profile.M-1)/2);

%% Compute cp
cp = 1-field.gamma(1:end-1).^2;

%% Plot profile
figure(); hold on; box on;
plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
quiver( profile.panels.centre.X,    profile.panels.centre.Y, ...
        -cp.*profile.panels.normal(:,1), -cp.*profile.panels.normal(:,2))
axis equal; xlabel('x'); ylabel('y')

% Plot centre of mass and pressure
centroid = computeCentroid(profile);
centreOfPressure = computeCentreOfPressure(profile,cp);

h1 = plot(centroid.x, centroid.y, 'ko', 'MarkerFaceColor', 'k');
h2 = plot(centreOfPressure.x, centreOfPressure.y, 'ro', 'MarkerFaceColor', 'r');

legend([h1, h2], ['Centre of mass = (', num2str(centroid.x,'%.2f'),', ' num2str(centroid.y,'%.2f'), ')'],...
                 ['Centre of pressure = (', num2str(centreOfPressure.x,'%.2f'),', ' num2str(centreOfPressure.y,'%.2f') ')'])
