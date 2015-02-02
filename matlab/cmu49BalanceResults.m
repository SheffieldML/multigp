function [err, errR] = cmu49BalanceResults(dataSetName, experimentNo)

% CMU49BALANCERESULTS descripton.
  
% MULTIGP
  
fillColor = [0.7 0.7 0.7];

load demCmu49BalanceArm5
model2=model;

global printDiagram

capName = dataSetName;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat'], 'model');


% load data
[y, void, yTest, void] = lvmLoadData(dataSetName);



% Get the time index.
fps = 120/32;
scaleVal = sqrt(sum(var(y)));
y = y/scaleVal;
X{1}{1} = [0];
Y{1}{1} = [0];
X{2}{1} = [0];
Y{2}{1} = [0];
for i = 1:size(y, 2)
  Y{1}{i+1} = y(1:35, i);
  Y{2}{i+1} = y(36:end, i);
  X{1}{i+1} = (1:35)'/fps;
  X{2}{i+1} = (1:30)'/fps;
end


% Get the time index.
fps = 120/32;
yTest = yTest/scaleVal;
Xtest{1} = [0];
Ytest{1} = [0];
for i = 1:3
  Ytest{i+1} = yTest(1:29, i);
  Xtest{i+1} = (1:29)'/fps;
end
for i = 4:9
  Ytest{i+1} = yTest(1, i);
  Xtest{i+1} = (1)'/fps;
end

xRTest = yTest(:, 1:3);

options = multigpOptions('ftc');
options.optimiser = 'conjgrad';
options.kernType = 'lfm';
options.nlf = 2;
options.tieOptions.selectMethod = 'free';

q = 1;
d = size(y, 2)+1;


testt = linspace(0, 10, 200)';
% for compNo = 1:length(model.comp)
%   [mu, varsig] = multigpPosteriorMeanVar(model.comp{compNo}, testt);
  
%   endVal = 0;
%   for i = 1:length(mu)
%     startVal = endVal + 1;
%     endVal = endVal + length(model.comp{compNo}.X{i});
%     figure, 
%     fill([testt; testt(end:-1:1)], ...
%          [mu{i}; mu{i}(end:-1:1)] ...
%          + 2*[sqrt(varsig{i}); -sqrt(varsig{i}(end:-1:1))], ...
%          fillColor,'EdgeColor',fillColor)
%     hold on
%     a = plot(testt, mu{i});
%     set(a, 'linewidth', 2);
  
%     b = plot(model.comp{compNo}.X{i}, model.comp{compNo}.y(startVal:endVal, ...
%                                                       1), '.');
%     set(b, 'markersize', 20);
%   end
% end



fps = 120/32;
ttest = (1:29)'*fps;

testModel = multigpCreate(q, d, Xtest, Ytest, options);
param = modelExtractParam(model);
testModel = modelExpandParam(testModel, param);
testModel.mu =  [0 yTest(1, :)];
testModel.B = yTest(1, :).*testModel.mu_D;
param = modelExtractParam(testModel);
testModel = modelExpandParam(testModel, param);
testt = linspace(0, 9, 100)';
[mu, varsig] = multigpPosteriorMeanVar(testModel, testt);
predVal = multigpPosteriorMeanVar(testModel, Xtest{2});

endVal = 0;
figNo = 1;
for i = 1:length(mu)
  startVal = endVal + 1;
  endVal = endVal + length(testModel.X{i});
  figure(figNo);
  fill([testt; testt(end:-1:1)], ...
       [mu{i}; mu{i}(end:-1:1)]*scaleVal ...
       + 2*[sqrt(varsig{i}); -sqrt(varsig{i}(end:-1:1))]*scaleVal, ...
       fillColor,'EdgeColor',fillColor)
  hold on
  a = plot(testt, mu{i}*scaleVal);
  set(a, 'linewidth', 2);
  if i>1
    b = plot(Xtest{2}, yTest(:, i-1)*scaleVal, '.');
    set(b, 'markersize', 20);
    if i > 4
      [muR, varsigmaR] = gpPosteriorMeanVar(model2{i-4}, xRTest);
      c = errorbar(Xtest{2}, muR*scaleVal, sqrt(varsigmaR)*scaleVal, 'x');
      set(c, 'linewidth', 2)
      errR(i-4) = sqrt(mean((muR - yTest(:, i-1)).^2))*scaleVal;
      err(i-4) = sqrt(mean((predVal{i} - yTest(:, i-1)).^2))*scaleVal;
    end
  end
  %  set(gca, 'xtick', [-2 -1 0 1 2]);
  %  set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
  ylim = get(gca, 'ylim');
  ylim(1) = min([floor(ylim(1)) 0]);
  ylim(2) = max([ceil(ylim(2)) 0]);
  set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [0 9], 'ylim', ylim)
  if exist('printDiagram') & printDiagram
    printPlot(['dem' capName num2str(experimentNo) '_' num2str(figNo)], '../tex/diagrams', '../html');
  end
  figNo = figNo + 1;

end

  
end
