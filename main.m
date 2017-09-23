addpath('./panel/')
addpath('./geometry/')

% Parameters
profile.c = 1;
profile.N = 100;
NACA = [4 4 1 2];

% Create profile and panels
profile = naca4(profile,[4,4,1,2]);
profile = create_panels(profile);

% Plot profile
figure(); hold on; box on;
plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
axis equal; xlabel('x'); ylabel('y')

% Solve potential flow
field=potential(profile);