function h = fxDataHintonDiag(dataSetName, kernType, approx, experimentNo, varargin);

% FXDATAHINTONDIAG Show the Hinton diagram of the results of any given
% demo for the foreign exchange rates data set.
%
% FORMAT
% DESC Shows the results of the given demo for the FX data set.
% RETURN
% ARG dataSetName : The data set to load.
% ARG kernType : The kernel used in the simulation ('sim' or 'simwhite').
% ARG approx : The approximation used ('ftc' or 'dtc').
% ARG experimentNo : Experiment number used in the demo.
%
% FORMAT Does the same as above but allows to specify a base directory for
% loading the result's file and flags regarding the saving and display of
% the results.
% RETURN
% ARG dataSetName : The data set to load.
% ARG kernType : The kernel used in the simulation ('sim' or 'simwhite').
% ARG approx : The approximation used ('ftc' or 'dtc').
% ARG experimentNo : Experiment number used in the demo.
% ARG baseDirResults : Base directory for the result's file.
% ARG flag : Binary variable specifying whether the figure should be saved.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, demSimwhiteDtcFxData, loadFxData

% MULTIGP


if (nargin >= 5) & ~isempty(varargin{1})
    baseDirResults = varargin{1};
else
    baseDirResults = './';
end
if (nargin == 6)
    saveFigure = varargin{end};
else
    saveFigure = 0;
end

offsetX = 1;

capName = dataSetName;
capName(1) = upper(capName(1));
kernName = kernType;
kernName(1) = upper(kernType(1));
approxName = approx;
approxName(1) = upper(approxName(1));

load([baseDirResults 'dem' kernName approxName capName 'Exp' num2str(experimentNo) '.mat'], ...
    'model', 'missingData', 'meanVal', 'scaleVal');

nlf = model.comp{1}.nlf;
nout = model.comp{1}.nout;
[params, names] = modelExtractParam(model);
S = zeros(nout, nlf);
for i=1:nlf
    if strcmp(model.comp{1}.kern.comp{i}.comp{i}.type, 'rbf')
        for j=1:nout
            S(j, i) = model.comp{1}.kern.comp{i}.comp{nlf+j}.variance;
        end
    else
        for j=1:nout
            S(j, i) = model.comp{1}.kern.comp{i}.comp{nlf+j}.sensitivity;
        end
    end
end

scaleFactor = std(S);
Snorm = S./repmat(scaleFactor, nout, 1);
permutation = [1 2 3 5 7 8 4 13 6 10 12 9 11];
Spermute = Snorm(permutation, :);
[xvals, yvals, color] = hintmat(Spermute');
h = hinton(Spermute');
xTick = mean(xvals(1:nlf:end,:).');
yTick = mean(yvals(1:nlf,:).');
xLim = [xTick(1)-0.5 xTick(end)+0.5];
yLim = [yTick(1)-0.5 yTick(end)+0.5];
set(gca, 'Visible', 'on', 'Box', 'on', 'Color', [0.5 0.5 0.5]);
set(gca, 'Units', 'normalized', 'ActivePositionProperty', 'OuterPosition', 'OuterPosition', [0 0 1 1]);
set(gca, 'XLim', xLim, 'YLim', yLim);
set(gca, 'XTick', xTick, 'YTick', yTick);
set(gca, 'XTickLabel', ['XAU'; 'XAG'; 'XPT'; 'EUR'; 'GBP'; 'CHF'; 'CAD'; 'MXN'; 'JPY'; 'HKD'; 'KRW'; 'AUD'; 'NZD']);
set(gca, 'YTickLabel', ['LF1'; 'LF2'; 'LF3'; 'LF4']);

if saveFigure
    set(h, 'InvertHardCopy', 'off')
    fileName = ['dem' kernName '_' capName '_Exp' num2str(experimentNo) '_HintonPlot'];
    print('-depsc', [baseDirResults fileName]);
    saveas(gcf,[baseDirResults fileName], 'fig');
    pos = get(gcf, 'paperposition');
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth);
    print('-dpng', [baseDirResults fileName])
    set(gca, 'lineWidth', lineWidth);
    set(gcf, 'paperposition', origpos);
end
