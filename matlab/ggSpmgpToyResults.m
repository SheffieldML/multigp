function ggSpmgpToyResults(dataSetName, experimentNo, XTemp, ...
    yTemp, XGT, fGT, initialLoc)

% GGSPMGPTOYRESULTS Show the prediction results for the demSpmgpGgToy's demos.
% FORMAT
% DESC shows the prediction results for the demSpmgpGgToy's
% demos.
% ARG dataSetName : name of the dataset used for the demo
% ARG experimentNo : number of the experiment
% ARG XTemp : input locations training data
% ARG yTemp : output values for training data
%
% DESC shows the prediction results for the demSpmgpGgToy's
% demos.
% ARG dataSetName : name of the dataset used for the demo
% ARG experimentNo : number of the experiment
% ARG XTemp : input locations training data
% ARG yTemp : output values for training data
% ARG initialLoc : initial location of spseudo inputs.
%
% COPYRIGHT : Mauricio Alvarez, 2008

% MULTIGP

capName = dataSetName;
capName(1) = upper(capName(1));
load(['demSpmgp' capName num2str(experimentNo) '.mat'], 'model');

fontsize = 25;
linewidth = 3;
markersize = 20;

saveFigures = false;
scaleVal = 1;
X = cell(size(yTemp, 2),1);
y = cell(size(yTemp, 2),1);

for i = 1:size(yTemp, 2)
    y{i} = yTemp{i};
    X{i} = XTemp{i};
end

Xt = linspace(min(X{model.nlf+1})-0.5,max(X{model.nlf+1})+0.5,200)';
%Xt = linspace(-1.3,1.3,200)';
%Xt = [-1.25 1.27];
[mu, varsigma] = multigpPosteriorMeanVar(model, Xt);
close all
%xlim = [min(model.X_u)  max(model.X_u)];
%xlim = [-1.2 1.2];
%ylim = [-12 12];
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
    if strcmp(model.kernType, 'lmc')
        muP = mu{k};
        varsigmaP = varsigma{k};
        f = [(muP+2*real(sqrt(varsigmaP)))*scaleVal;flipdim((muP-2*real(sqrt(varsigmaP)))*scaleVal,1)];
        a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
        a =[ a plot(Xt, muP*scaleVal,'k-')];
        if k>model.nlf*model.rankCorregMatrix
            c =plot(X{k-model.nlf*model.rankCorregMatrix},y{k-model.nlf*model.rankCorregMatrix}*scaleVal,'k.');
            d =plot(XGT{k-model.nlf*model.rankCorregMatrix}, fGT{k-model.nlf*model.rankCorregMatrix}, 'k--');
        end
        minimum = min((muP-2*real(sqrt(varsigmaP)))*scaleVal);
        maximum = max((muP+2*real(sqrt(varsigmaP)))*scaleVal);
%         if isfield(model, 'X_u') && ~isempty(model.X_u);
%             b = plot(model.X_u, -10, 'kx');
%             set(b, 'linewidth', linewidth)
%             set(b, 'markersize', markersize);
%             if nargin ==7
%                 hold on
%                 d = plot(initialLoc, maximum*0.9, 'kx');
%                 set(d, 'linewidth', linewidth)
%                 set(d, 'markersize', markersize);
%             end
%         end
        if k>model.nlf*model.rankCorregMatrix
            set(c,   'markersize', 0.7*markersize);
            set(d,   'lineWidth', linewidth);
        else
        end
    else
        f = [(mu{k}+2*real(sqrt(varsigma{k})))*scaleVal;flipdim((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal,1)];
        a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
        a =[ a plot(Xt, mu{k}*scaleVal,'k-')];
        if k>model.nlf
            c =plot(X{k-model.nlf},y{k-model.nlf}*scaleVal,'k.');
            d =plot(XGT{k-model.nlf}, fGT{k-model.nlf}, 'k--');
        end

        minimum = min((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal);
        maximum = max((mu{k}+2*real(sqrt(varsigma{k})))*scaleVal);
        if isfield(model, 'X_u') && ~isempty(model.X_u);
            b = plot(model.X_u, -10, 'kx');
            set(b, 'linewidth', linewidth)
            set(b, 'markersize', markersize);
            if nargin ==7
                hold on
                d = plot(initialLoc, maximum*0.9, 'kx');
                set(d, 'linewidth', linewidth)
                set(d, 'markersize', markersize);
            end
        end
        if k>model.nlf
            set(c,   'markersize', 0.7*markersize);
            set(d,   'lineWidth', linewidth);
        else
        end
    end
    %     ylabel('y', 'fontsize',fontsize);
    %     prop = get(g);
    %     poslabel = prop.Position;
    %     poslabel(1) = poslabel(1) -0.001*poslabel(1);
    %     ylabel('PEV', 'fontsize',fontsize, 'position', poslabel);
    %     g = xlabel('Input', 'fontsize',fontsize);
    %     prop = get(g);
    %     poslabel = prop.Position;
    %     poslabel(1) = 0;
    %     poslabel(2) = -14;
    %     xlabel('Input', 'fontsize',fontsize, 'position', poslabel);
    set(a,   'lineWidth', 2);
    %set(gca, 'fontname', 'arial', 'fontsize', fontsize, 'xlim', xlim, 'ylim', ylim, 'Color', 'none')
    set(gca, 'Color', 'none')
    %set(gca, 'position', [0.06 0.08 0.9 0.9])
    %high = get(0, 'screensize');
    %set(gcf, 'position', high)
    %set(gcf, 'PaperPositionMode', 'auto');
    box on
    if saveFigures
        fileName = ['toy1D' upper(model.approx) num2str(k-model.nlf)];
        print('-dpdf', ['./resultsToy1D/' fileName]);
        print('-depsc', ['./resultsToy1D/' fileName], '-loose');
        print('-dpng', ['./resultsToy1D/' fileName]);
    end
end
