% FXDATACOMPAREVARIATIONALBOUND description

% MULTIGP

% Data

R = 1:1:6;
Fv_best = [-2115.52 -1173.91 -692.28 -235.11 216.31 555.25];
Fv_sim = [-2115.52 -1173.91 -692.28 -257.39 216.31 555.25];
Fv_ou = [-2257.29 -1408.13 -977.31 -440.58 -123.49 95.83];

baseDirResults = './';

% Plotting the results

figure
hold on
% Plotting the true data
c = plot(R, Fv_best, 'k-*');
set(c, 'markerSize', 10, 'lineWidth', 2);
d = plot(R, Fv_sim, 'r-s');
set(d, 'markerSize', 10, 'lineWidth', 2);
e = plot(R, Fv_ou, 'b-o');
set(e, 'markerSize', 10, 'lineWidth', 2);
set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
hl = legend('Best', 'Only SIM', 'Only OU', 'Location', 'SouthEast');
box on
fileName = 'fxDataCompareVariationalBound';
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
