close all
clear all

addpath('./panel/')
addpath('./geometry/')
addpath('./utilities/')
set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');


% Parameters
profile.c = 1;
profile.N = 40; %(Oben und unten N-1 Panels + Hinterseite) -> M=1+2(N-1)=2N-1 Randelemente
profile.alfa = 2*pi/180;
NACA = [4 4 1 2];

% Create profile and panels
profile = naca4(profile,[4,4,1,2]);
profile = create_panels(profile);

% Plot profile
% figure(); hold on; box on;
% plot([profile.panels.X]',[profile.panels.Y]','k','Linewidth',2);
% axis equal; xlabel('x'); ylabel('y')

% Solve potential flow
field=potential(profile); Nle = round((profile.M-1)/2);
figure()
hold on; box on
plot(mean(profile.panels.X(1:Nle,:),2), (1-field.gamma(1:Nle).^2));
plot(mean(profile.panels.X(Nle:end,:),2), (1-field.gamma(Nle:end-1).^2));
load naca4412_alfa2.txt
plot(naca4412_alfa2(:,1),naca4412_alfa2(:,3),'.');
ylim([-1.2 1]);
xlim([0 1]);
legend('$C_p=1-\gamma(s)^2$ Saugseite','$C_p=1-\gamma(s)^2$ Druckseite','Vergleichswerte');

