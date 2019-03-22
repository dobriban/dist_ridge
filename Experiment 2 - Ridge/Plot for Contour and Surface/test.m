%plot are
cd('C:\Dropbox\Projects\Distributed Ridge\Experiments\Experiment 2 - Ridge\Plots');
addpath '..\..\'
%%

%% Plot ARE
alpha = linspace(0,5,30);
gamma = linspace(0,3,30);
[X,Y] = meshgrid(alpha, gamma);
Z = are(X,Y);

figure
colormap(brighten(jet,0.6))
surf(X,Y,Z)
xlim([min(alpha),max(alpha)])
ylim([min(gamma),max(gamma)])
zlim([0,max(max(Z))])
xlabel('\alpha');
ylabel('\gamma');
zlabel('ARE');
%%
saveTightFigure(gcf,'ARE_Surface.pdf')

%% Contour

figure
contour(X,Y,Z,'ShowText','on')
xlim([min(alpha),max(alpha)])
ylim([min(gamma),max(gamma)])
xlabel('\alpha');
ylabel('\gamma');
colormap(brighten(jet,0.6))
%%
saveTightFigure(gcf,'ARE_Contour.pdf')

