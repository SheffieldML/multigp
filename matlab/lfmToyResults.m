function lfmToyResults(dataSetName, experimentNo)

% LFMTOYRESULTS description.
  
% MULTIGP
  
capName = dataSetName;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat'], 'model');
saveFigures = false;
% load data
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Set equal variances 
scaleVal = sqrt(sum(var(yTemp)));
yTemp = yTemp/scaleVal;
yTestTemp = yTestTemp/scaleVal;

X = cell(size(yTemp, 2)+ model.nlf,1);
y = cell(size(yTemp, 2)+ model.nlf,1);
XTest = cell(size(yTemp, 2)+ model.nlf,1);
yTest = cell(size(yTemp, 2)+ model.nlf,1);

for i = 1:model.nlf
    X{i} = 0;
    y{i} = 0;
end
for i = 1:size(yTemp, 2)
  y{i+model.nlf} = yTemp(:, i);
  X{i+model.nlf} = XTemp;
end

for i = 1:model.nlf,
    XTest{i} = 0;
    yTest{i} = 0;
end
for i = 1:size(yTemp, 2)
    yTest{i+model.nlf} = yTestTemp(:, i);
    XTest{i+model.nlf} = XTestTemp;
end


Xt = linspace(min(X{model.nlf+1})-0.2,max(X{model.nlf+1})+0.2,100)';
[mu, varsigma] = multigpPosteriorMeanVar(model, Xt);
close all

if strcmp(model.approx,'ftc')
    nFigs = model.d;
    xlim = [min(X{model.nlf+1}) max(X{model.nlf+1})];
else
    nFigs = model.nout+model.nlf;
    xlim = [min(model.X_u)  max(model.X_u)];
end
for k=1:nFigs,
    figure
    hold on
    f = [(mu{k}+2*real(sqrt(varsigma{k})))*scaleVal;flipdim((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal,1)];
    a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(Xt, mu{k}*scaleVal,'k-')];
    c =plot(X{k},y{k}*scaleVal,'k.');
    minimum = min((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal);
    maximum = max((mu{k}+2*real(sqrt(varsigma{k})))*scaleVal);
    if isfield(model, 'X_u') && ~isempty(model.X_u);
        b = plot(model.X_u, minimum*0.9, 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
%    ylim = [min(minimum,min(y{k}*scaleVal)) max(maximum,max(y{k}*scaleVal))];
    set(a,   'lineWidth', 2);
    set(c,   'markersize', 10);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim,  'Color', 'none')
    box on
    if saveFigures==1
        fileName = ['Toy_prediction' 'full' num2str(k)];
        print('-depsc', ['./results/' fileName]);
        saveas(gcf,['./results/' fileName],'fig');
        pos = get(gcf, 'paperposition');
        origpos = pos;
        pos(3) = pos(3)/2;
        pos(4) = pos(4)/2;
        set(gcf, 'paperposition', pos);
        lineWidth = get(gca, 'lineWidth');
        set(gca, 'lineWidth', lineWidth);
        print('-dpng', ['./results/' fileName])
        set(gca, 'lineWidth', lineWidth);
        set(gcf, 'paperposition', origpos);
    end
end
%err(i-4) = sqrt(mean((predVal{i} - yTest(:, i-1)).^2))*scaleVal;

