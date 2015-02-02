function ggToyResults(dataSetName, experimentNo, XTemp, yTemp, XGT, fGT)

% GGTOYRESULTS Show the prediction results for the demGgToy demo.
% FORMAT
% DESC Show the prediction results for the demGgToy demo.
% ARG dataSetName : name of the dataset used for the demo
% ARG experimentNo : number of the experiment
% ARG XTemp : input locations training data
% ARG yTemp : output values for training data
%
% COPYRIGHT : Mauricio Alvarez, 2008

% MULTIGP


capName = dataSetName;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat'], 'model');
saveFigures = false;
scaleVal = 1;

fontsize = 25;
linewidth = 3;
markersize = 20;

X = cell(size(yTemp, 2)+model.nlf,1);
y = cell(size(yTemp, 2)+model.nlf,1);

for j=1:model.nlf
    y{j} = 0;
    X{j} = 0;
end
for i = 1:size(yTemp, 2)
    y{i+model.nlf} = yTemp{i};
    X{i+model.nlf} = XTemp{i};
end

Xt = linspace(min(X{model.nlf+1})-0.5,max(X{model.nlf+1})+0.5,200)';
Xt = linspace(-1.3,1.3,200)';
%Xt = [-1.25 1.27];
[mu, varsigma] = multigpPosteriorMeanVar(model, Xt);
% close all
% xlim = [min(model.X_u)  max(model.X_u)];
% xlim = [-1.2 1.2];
% ylim = [-12 12];

close all
%
%xlim = [-1.05 1.05];
if strcmp(model.kernType, 'lmc')
    nFigs = model.nout+model.nlf*model.rankCorregMatrix;
    muP = cell(nFigs,1);
    varsigmaP = cell(nFigs,1);
    for r=1:model.nlf
        for j=1:model.rankCorregMatrix
            linI = (r-1)*model.rankCorregMatrix + j;
            muP{linI} = mu{r}{j};
            varsigmaP{linI} = varsigma{r}{j};
        end
    end
    muP(model.nlf*model.rankCorregMatrix+1:end) = mu(model.nlf+1:end);
    varsigmaP(model.nlf*model.rankCorregMatrix+1:end) = varsigma(model.nlf+1:end);
    mu = muP;
    varsigma = varsigmaP;
    X = XTemp;
    y = yTemp;
else
    if strcmp(model.approx,'ftc')
        nFigs = model.d;
    else
        nFigs = model.nout+model.nlf;
    end
end
for k=1:nFigs,
    figure
    hold on
    f = [(mu{k}+2*real(sqrt(varsigma{k})))*scaleVal;flipdim((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal,1)];
    a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(Xt, mu{k}*scaleVal,'k-')];
    if strcmp(model.kernType, 'lmc')
        if k>(model.nlf*model.rankCorregMatrix)
            c =plot(X{k-model.nlf*model.rankCorregMatrix},y{k-model.nlf*model.rankCorregMatrix}*scaleVal,'k.');
            d =plot(XGT{k-model.nlf*model.rankCorregMatrix}, fGT{k-model.nlf*model.rankCorregMatrix}, 'k--');
        end
        minimum = min((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal);
        maximum = max((mu{k}+2*real(sqrt(varsigma{k})))*scaleVal);
        if isfield(model, 'X_u') && ~isempty(model.X_u);
            b = plot(model.X_u, -8, 'kx');
            set(b, 'linewidth', 2)
            set(b, 'markersize', 10);
        end
        if k>model.nlf*model.rankCorregMatrix
            %ylim = [min(minimum,min(y{k}*scaleVal)) max(maximum,max(y{k}*scaleVal))];
            set(c,   'markersize', 0.7*markersize);
            set(d,   'lineWidth', linewidth);
        else
            %ylim = [min(minimum,min(y{model.nlf+1}*scaleVal)) max(maximum,max(y{model.nlf+1}*scaleVal))];
        end
    else
        if k>model.nlf
            c =plot(X{k},y{k}*scaleVal,'k.');
            d =plot(XGT{k-model.nlf}, fGT{k-model.nlf}, 'k--');
        end
        minimum = min((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal);
        maximum = max((mu{k}+2*real(sqrt(varsigma{k})))*scaleVal);
        if isfield(model, 'X_u') && ~isempty(model.X_u);
            b = plot(model.X_u, -8, 'kx');
            set(b, 'linewidth', 2)
            set(b, 'markersize', 10);
        end
        if k>model.nlf
            %ylim = [min(minimum,min(y{k}*scaleVal)) max(maximum,max(y{k}*scaleVal))];
            set(c,   'markersize', 0.7*markersize);
            set(d,   'lineWidth', linewidth);
        else
            %ylim = [min(minimum,min(y{model.nlf+1}*scaleVal)) max(maximum,max(y{model.nlf+1}*scaleVal))];
        end
    end

    %        ylabel('y', 'fontsize',fontsize);
    %     prop = get(g);
    %     poslabel = prop.Position;
    %     poslabel(1) = poslabel(1) -0.001*poslabel(1);
    %     ylabel('PEV', 'fontsize',fontsize, 'position', poslabel);
    %    g = xlabel('Input', 'fontsize',fontsize);
    %    prop = get(g);
    %    poslabel = prop.Position;
    %    poslabel(1) = 0;
    %    poslabel(2) = -14;
    %    xlabel('Input', 'fontsize',fontsize, 'position', poslabel);
    set(a,   'lineWidth', 2);
    set(gca, 'fontname', 'arial', 'fontsize', fontsize, 'xlim', xlim, 'ylim', ylim, 'Color', 'none')
    %set(gca, 'position', [0.06 0.08 0.9 0.9])
    %    high = get(0, 'screensize');
    %    set(gcf, 'position', high)
    %    set(gcf, 'PaperPositionMode', 'auto');
    box on
    if saveFigures
        if strcmp(model.kernType, 'lmc')
            fileName = ['toy1D' upper(model.approx) num2str(k-model.nlf*model.rankCorregMatrix)];
        else
            fileName = ['toy1D' upper(model.approx) num2str(k-model.nlf)];
        end
        print('-dpdf', ['./resultsToy1D/' fileName]);
        print('-depsc', ['./resultsToy1D/' fileName], '-loose');
        print('-dpng', ['./resultsToy1D/' fileName]);
    end
end
