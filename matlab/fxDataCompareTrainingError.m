% FXDATACOMPARETRAININGERROR description

% MULTIGP

% Data

R = 1:1:6;
smse_train_best = [33.15 19.66 11.61 6.49 5.34 4.50];
smse_train_sim = [34.88 20.29 13.97 11.15 7.16 5.04];
smse_train_ou = [33.15 19.66 13.28 7.94 5.34 4.50];
smse_train_lmc = [32.71 19.66 13.32 7.46 5.18 3.83 2.67];

baseDirResults = './';

% Plotting the results

figure
hold on
% Plotting the true data
c = plot(R, smse_train_best, 'k-*');
set(c, 'markerSize', 10, 'lineWidth', 2);
d = plot(R, smse_train_sim, 'r-s');
set(d, 'markerSize', 10, 'lineWidth', 2);
e = plot(R, smse_train_ou, 'b-o');
set(e, 'markerSize', 10, 'lineWidth', 2);
f = plot(R, smse_train_lmc(1:length(R)), 'g-v');
set(f, 'markerSize', 10, 'lineWidth', 2);
set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
hl = legend('Best', 'Only SIM', 'Only OU', 'LMC', 'Location', 'NorthEast');
box on
fileName = 'fxDataCompareTrainingError';
print('-depsc', [baseDirResults 'results/' fileName]);
saveas(gcf,[baseDirResults 'results/' fileName],'fig');
pos = get(gcf, 'paperposition');
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth);
print('-dpng', [baseDirResults 'results/' fileName])
set(gca, 'lineWidth', lineWidth);
set(gcf, 'paperposition', origpos);
