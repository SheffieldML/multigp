% Script used to plot Figures 12 and 13 of the "Updated Results for the
% Foreign Exchange Data" technical report (NIPS 2009), which compare the
% performance of the different models considered.
%
% COPYRIGHT : David Luengo, 2009

% MULTIGP

close all
clc

R = 1:6;

FvMax = [-2117.4625 -1100.7633 -648.6683 -156.3234 227.9023 565.3821];
FvSmooth = [-2117.4625 -1175.5964 -711.0210 -258.1882 223.9023 565.3821];
FvNonSmooth = [-2251.8924 -1406.1206 -975.2001 -443.4104 130.8854 95.8532];

errorTrainMin = [33.1054 17.9501 11.7775 7.6471 5.2183 4.1128];
errorTrainBestBound = [34.8711 17.9501 11.7775 7.6471 7.0662 5.1186];
errorTrainSmooth = [34.8711 20.3079 14.0653 11.1231 7.4381 5.1186];
errorTrainNonSmooth = [33.1054 19.5480 13.2803 7.9073 5.3212 4.5062];
errorTrainPpca = [32.7148 19.6622 13.3207 7.4639 5.1798 3.8279];

errorTestMin = [53.4085 28.7506 31.2943 27.7547 19.3322 25.8079];
errorTestBestBound = [53.4085 28.7506 31.2943 32.2393 33.6536 100.0306];
errorTestSmooth = [53.4085 34.0249 33.2135 27.7547 37.9714 100.0306];
errorTestNonSmooth = [55.9586 30.3053 33.2843 34.4682 20.28660 70.3651];
errorTestPpca = [56.4046 39.2658 43.1681 43.1992 48.1880 56.4057];

figure(1)
h = plot(R, FvSmooth, 'rs-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
hold on
h = plot(R, FvNonSmooth, 'bo-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, FvMax, 'k*-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
xlabel('R', 'FontSize', 15, 'FontName', 'arial')
ylabel('Variational Bound', 'FontSize', 15, 'FontName', 'arial')
h = legend('Only SIM', 'Only OU', 'Best combination', 'Location', 'SouthEast');
set(h, 'FontSize', 15, 'FontName', 'arial', 'Color', 'none')
set(gca, 'XTick', R, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none');
box on
set(1, 'Position', [-3 35 1280 696])
print -depsc FxDataCompareVariationalBound.eps

figure(2)
hold on
h = plot(R, errorTrainSmooth, 'rs-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTrainNonSmooth, 'bo-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTrainPpca, 'mx-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTrainBestBound, 'gv-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTrainMin, 'k*-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
xlabel('R', 'FontSize', 15, 'FontName', 'arial')
ylabel('SMSE (Training)', 'FontSize', 15, 'FontName', 'arial')
h = legend('Only SIM', 'Only OU', 'PPCA', 'Best bound', 'Minimum error', 'Location', 'NorthEast');
set(h, 'FontSize', 15, 'FontName', 'arial', 'Color', 'none')
set(gca, 'XTick', R, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none');
box on
set(2, 'Position', [-3 35 1280 696])
print -depsc FxDataCompareTrainingError.eps

figure(3)
hold on
h = plot(R, errorTestSmooth, 'rs-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTestNonSmooth, 'bo-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTestPpca, 'mx-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTestBestBound, 'gv-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
h = plot(R, errorTestMin, 'k*-');
set(h, 'markerSize', 10);
set(h, 'lineWidth', 2)
xlabel('R', 'FontSize', 15, 'FontName', 'arial')
ylabel('SMSE (Test)', 'FontSize', 15, 'FontName', 'arial')
h = legend('Only SIM', 'Only OU', 'PPCA', 'Best bound', 'Minimum error', 'Location', 'NorthWest');
set(h, 'FontSize', 15, 'FontName', 'arial', 'Color', 'none')
set(gca, 'XTick', R, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none');
box on
set(3, 'Position', [-3 35 1280 696])
print -depsc FxDataCompareTestError.eps
