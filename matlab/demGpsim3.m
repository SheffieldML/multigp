% DEMGPSIM3 Demo of full multi output GP with the sim kernel using PITC
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process with the sim kernel using
%      PITC

% MULTIGP

clc
clear
rand('seed',1e6);
randn('seed',1e6);
%
dataSetName = 'dataBarencoOption_1_Genes_9';
experimentNo = 1;
ntrainX =200;
ntrainX2 =5;
iters =500;
approx = 'pitc';
saveFigures=1;
% load data
data = multigpLoadData(dataSetName);
data.nin = 2;  % Number of latent functions
missingData = cell(data.nout,1);
%
options.isSparse = 1;      % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;      % Indicates if the inputs are to be standarized (1) or not(0)
options.stdBiasY = 0;      % Indicates if the outputs are to be standarized (1) or not(0)
options.outKernelName = 'sim'; % Indicates the name of the output kernel: lfm, sim or gg

% Setup model
[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsSim(options, data, ntrainX, missingData, ntrainX2, approx);
options.optimiser = 'scg';
% Creates the model
model = multigpCreate(MYtrain, MXtrain, options);

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
end

Zt = linspace(min(Xtrain{1})-0.2,max(Xtrain{1})+0.2,100)';
% This part is to plot the posterior distribution
[mean_u, varsigma_u] = multigpPosteriorMeanVar(model, Zt);
% xlim = [min(model.X_u) max(model.X_u)];
for j =1:model.nlf
    figure    
    hold on
    f = [mean_u{j}+2*real(sqrt(varsigma_u{j}));flipdim(mean_u{j}-2*real(sqrt(varsigma_u{j})),1)];
    a = fill([Zt; flipdim(Zt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
 %   a =[ a plot(data.Ztrain{j}, data.Utrain{j}, 'k--')];
    a =[ a plot(Zt,mean_u{j},'k-')];
    set(a,   'lineWidth', 2);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none')
    if saveFigures==1
        fileName = ['Toy_posterior' 'full' num2str(j)];
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
