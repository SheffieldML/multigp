% FXDATACOMPARETESTERROR description

% MULTIGP

% Data

R = 1:1:6;
smse_test_best = [53.61 29.98 33.00 27.95 20.33 28.47];
smse_test_sim = [53.61 34.06 35.48 28.15 41.92 110.17];
smse_test_ou = [56.00 29.98 33.00 34.89 20.33 70.37];
smse_test_lmc = [56.41 39.27 43.17 43.20 48.19 56.41 59.75];

baseDirResults = './';

% Plotting the results

figure
hold on
% Plotting the true data
c = plot(R, smse_test_best, 'k-*');
set(c, 'markerSize', 10, 'lineWidth', 2);
d = plot(R, smse_test_sim, 'r-s');
set(d, 'markerSize', 10, 'lineWidth', 2);
e = plot(R, smse_test_ou, 'b-o');
set(e, 'markerSize', 10, 'lineWidth', 2);
f = plot(R, smse_test_lmc(1:length(R)), 'g-v');
set(f, 'markerSize', 10, 'lineWidth', 2);
set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
hl = legend('Best', 'Only SIM', 'Only OU', 'LMC', 'Location', 'NorthWest');
box on
fileName = 'fxDataCompareTestError';
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
