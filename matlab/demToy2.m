% DEMTOY2 Demonstration of sparse multioutput model using FITC.

% In this demo, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.
% When changing the kernel, the fix values in this code and the indeces in
% the tieParam vector in multigpCreate must also be changed.


% MULTIGP

clc
clear
rand('seed',1e6);
randn('seed',1e6);
%
dataSetName = 'sample1D41';
experimentNo = 1;
ntrainX =200;
ntrainX2 =30;
iters =1000;
approx = 'fitc';
saveFigures=1;
% load data
data = multigpLoadData(dataSetName);
data.nin = 1;  % Number of latent functions
missingData = cell(data.nout,1);
missingData{1} = 101:160;
missingData{4} = 21:90;

options.isSparse = 1;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 0;  % Indicates if the outputs are to be standarized
options.outKernelName = 'gg'; % Indicates the name of the output kernel: lfm, sim or gg
% Setup model
[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsToy(options, data, ntrainX,missingData, ntrainX2, approx);
% Specify the options for the kernel and tie the parameters
options.optimiser = 'scg';
% Creates the model
model = multigpCreate(MYtrain, MXtrain, options);

% Trains the model and counts the training time
ini_t = cputime;
[model, params] = multigpOptimise(model,1,iters);
total_t = cputime - ini_t;

% This part is to show the plots of the mean prediction and error bars
Xt = linspace(min(Xtrain{1})-0.2,max(Xtrain{1})+0.2,500)';
[mu, varsigma] = multigpPredictionMeanVar(model, Xt, options);
close all
xlim = [-1.2 1.2];
for k=1:model.nout,
    figure
    hold on
    f = [mu{k}+2*real(sqrt(varsigma{k}));flipdim(mu{k}-2*real(sqrt(varsigma{k})),1)];
    a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(data.x, data.F(:,k), 'k--')];
    a =[ a plot(Xt, mu{k},'k-')];
    c =plot(Xtrain{k},Ytrain{k},'k.');
    minimum = min(mu{k}-2*real(sqrt(varsigma{k})));
    maximum = max(mu{k}+2*real(sqrt(varsigma{k})));
    if isfield(model, 'X_u') && ~isempty(model.X_u);
        b = plot(model.X_u, -8, 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
    ylim = [-10 10];
    set(a,   'lineWidth', 2);
    set(c,   'markersize', 10);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'ylim', ylim, 'Color', 'none')
    box on
    if saveFigures==1
        fileName = ['Toy_prediction' approx num2str(k)];
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


if model.isSparse
    % This part is to plot the posterior distribution
    [mean_u, varsigma_u] = multigpPosteriorMeanVar(model);
    % xlim = [min(model.X_u) max(model.X_u)];
    for j =1:model.nlf
        figure
        [X_u, index] = sort(model.X_u);
        d = plot(X_u,mean_u{j}(index),'b-');
        hold on
        d = [d plot(X_u,mean_u{j}(index) + 2*real(sqrt(varsigma_u{j}(index))),'b--')];
        d = [d plot(X_u,mean_u{j}(index) -  2*real(sqrt(varsigma_u{j}(index))),'b--')];
        set(d,   'lineWidth', 2);
        set(gca, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none')
    end
end


