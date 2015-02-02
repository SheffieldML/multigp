% DEMMOCAP1 Demo of full multi output GP with the lfm kernel.
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process with the lfm kernel.


% MULTIGP

clc

if exist('current_seed.mat', 'file')
    load current_seed;
    rand('twister', seed);
    randn('state', seed);
else
    seed = 1e6;
    rand('seed',1e6);
    randn('seed',1e6);
end

if (~exist('dataSetName', 'var'))
    dataSetName = 'datasetMOCAPSubject49Train18Test19Dec32SelChan7-Balance-RightArm';
end

experimentNo = 1;
ntrainX = 35;
iters = 1000;
saveFigures = 0;

% Load data

data = multigpLoadData(dataSetName);

% Set the number of latent functions

if exist('numberLatentForces', 'var')
    data.nin = numberLatentForces;
else
    data.nin = 1;
end

if (~exist('useInd', 'var'))
    useInd = 1;
end

missingData = cell(data.nout,1);

options.isSparse = 0;      % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;      % Indicates if the inputs are to be standarized (1) or not(0)
options.stdBiasY = 1;      % Indicates if the outputs are to be standarized (1) or not(0)
options.outKernelName = 'lfm'; % Indicates the name of the output kernel: lfm, sim or gg
options.tieSelectMethod = 'free'; % Indicates the method to tie the parameters

% Setup model

[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsMocap(options, data, ntrainX, cell(data.nout,1), [], 'none', useInd);
options.optimiser = 'conjgrad';

% Creates the model

model = multigpCreate(MYtrain, MXtrain, options);

% Initialization

Delta_t = data.Xtrain{1}(end)-data.Xtrain{1}(1);
delta_t = mean(diff(data.Xtrain{1}));
% invWidthlfm = 1./(((1:data.nin)*Delta_t/data.nin).^2);
invWidthlfm = 1./(((1:data.nin)*Delta_t/(data.nin+1)).^2);
paramsOrig = kernExtractParam(model.kern);
paramsMod = paramsOrig;
paramsMod(4) = log(invWidthlfm(1));
paramsMod(4*data.nout+2:data.nout+1:(2+data.nin)*data.nout+data.nin) = log(invWidthlfm(2:data.nin));
if options.includeInd
    paramsMod((data.nin+3)*data.nout+data.nin+1:2:(data.nin+5)*data.nout+data.nin-1) = -2*log(delta_t);
%     paramsMod((data.nin+3)*data.nout+data.nin+2:2:(data.nin+5)*data.nout+data.nin) = log(3);
end
model.kern = kernExpandParam(model.kern, paramsMod);

% Trains the model and counts the training time

ini_t = cputime;
[model, params] = multigpOptimise(model,1,iters);
total_t = cputime - ini_t;

% This part is to show the plots of the mean prediction and error bars

Xt = linspace(min(Xtrain{1})-0.2,max(Xtrain{1})+0.2,100)';
[mu, varsigma] = multigpPredictionMeanVar(model, Xt, options);
close all
xlim = [min(Xtrain{1}) max(Xtrain{1})];
for k=1:model.nout,
    figure
    hold on
    f = [mu{k}+2*real(sqrt(varsigma{k}));flipdim(mu{k}-2*real(sqrt(varsigma{k})),1)];
    a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(Xt, mu{k},'k-')];
    %    a =[ a plot(data.XtrainFull{k}, data.Ftrain{k}, 'k--')];
    c =plot(Xtrain{k},Ytrain{k},'k.');
    minimum = min(mu{k}-2*real(sqrt(varsigma{k})));
    maximum = max(mu{k}+2*real(sqrt(varsigma{k})));
    if isfield(model, 'X_u') && ~isempty(model.X_u);
        b = plot(model.X_u, -8, 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
    ylim = [min(minimum,min(Ytrain{k})) max(maximum,max(Ytrain{k}))];
    set(a,   'lineWidth', 2);
    set(c,   'markersize', 10);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'ylim', ylim, 'Color', 'none')
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
    pause(1)
end

% This part is to plot the posterior distribution

Zt = linspace(min(Xtrain{1})-5,max(Xtrain{1})+5,100)';
[mean_u, varsigma_u] = multigpPosteriorMeanVar(model, Zt, options);
% xlim = [min(model.X_u) max(model.X_u)];
for j =1:model.nlf
    figure
    hold on
    f = [mean_u{j}+2*real(sqrt(varsigma_u{j}));flipdim(mean_u{j}-2*real(sqrt(varsigma_u{j})),1)];
    a = fill([Zt; flipdim(Zt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    %    a =[ a plot(data.Ztrain{end-j+1}, data.Utrain{end-j+1}, 'k--')];
    %    a =[ a plot(data.Ztrain{j}, data.Utrain{j}, 'k--')];
    a =[ a plot(Zt,mean_u{j},'k-')];
    set(a,   'lineWidth', 2);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none')
    if saveFigures==1
        fileName = ['Toy_posterior' 'full' num2str(j) '_outputs_' num2str(options.nout)];
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
    pause(1)
end
