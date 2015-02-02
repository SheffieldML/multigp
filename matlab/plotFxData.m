function plotFxData(dataSetName, varargin);

% PLOTFXDATA Plot the FX data and save the figures.
%
% FORMAT
% DESC Plot the time series for the financial data experiments and save the
% resulting figures.
% ARG dataSetName : The data set to load.
%
% FORMAT
% DESC Does the same as above but allows to specify the directory where the
% data are to be saved and some flags related to the saving of the data.
% ARG dataSetName : The data set to load.
% ARG baseDirResults : Base directory for the result's file.
% ARG flags : Binary vector specifying: (1) whether the results should be
% plotted on the screen or not. (2) whether the figures should be saved.

% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, demSimwhiteDtcFxData, loadFxData,
% fxDataResultsPpca

% MULTIGP


if (nargin >= 2) & ~isempty(varargin{1})
    baseDirResults = varargin{1};
else
    baseDirResults = './';
end
if (nargin == 3)
    plotFigures = varargin{end}(1);
    saveFigures = varargin{end}(2);
else
    plotFigures = 1;
    saveFigures = 1;
end

data = loadFxData(dataSetName);

offsetX = 1;
nout = size(data.yTest, 2);

% Preparing the test data that are going to be plotted later

XTestPlot = cell(1, 1);
yTestPlot = cell(1, 1);
count = offsetX + (0:size(data.XTest, 1)-1)';
for i = 1:nout
    ind = find(data.yTest{i} ~= 0.0);
    yTestPlot{1}{i} = data.yTest{i}(ind);
    XTestPlot{1}{i} = count(ind);
end

% Plotting the data

close all
xlim = [min(XTestPlot{1}{end}) max(XTestPlot{1}{end})];
for k = 1:nout
    if plotFigures
        figure(k)
        c = plot(XTestPlot{1}{k}, yTestPlot{1}{k}, 'k.');
        set(c, 'markerSize', 20);
        set(c, 'lineWidth', 2)
        set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
        box on
        if saveFigures
            fileName = ['FxData_' dataSetName '_Fig' num2str(k)];
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
        end
    end
end

return