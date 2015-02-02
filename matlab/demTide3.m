% DEMTIDE3 Example of Sparse Multi Output Gaussian Process using PITC over the
% heigth tide dataset of the sensor network. 

% MULTIGP

% In this example, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.
% When changing the kernel, the fix values in this code and the indexes in
% the tieParam vector in multigpCreate must also be changed.

clc
clear
rand('seed',1e5);
randn('seed',1e5);
%
dataSetName = 'data_tide2';
experimentNo = 1;
ntrainX =1000;
ntrainX2 =100;
iters =50;
approx = 'pitc';
saveFigures=1;
% load data
data = multigpLoadData(dataSetName);
data.nin = 1;  % Number of latent functions
data.nout = 4;
missingData = cell(data.nout,1);
missingData{1} = 201:400 ;
missingData{4} =  501:700;

options.isSparse = 1;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 1;  % Indicates if the outputs are to be standarized
options.outKernelName = 'gg'; % Indicates the name of the output kernel: lfm, sim or gg

% Setup model
[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsTide(options, data, ntrainX,missingData, ntrainX2, approx );

options.optimiser = 'scg';

% Creates the model
model = multigpCreate(MYtrain, MXtrain, options);

% Trains the model and counts the training time
ini_t = cputime;
[model, params] = multigpOptimise(model,1,iters);
total_t = cputime - ini_t;

% This part is to show the plots of the mean prediction and error bars

[mu, varsigma] = multigpPredictionMeanVar(model, Xtest, options);
close all
xlim = [0 3];
sensors = {'Bramblemet', 'Chimet', 'Sotonmet', 'Cambermet'};
for k=1:model.nout,
    figure
    hold on
    [Xtestsort, indx] = sort(Xtest{k});
    f = [mu{k}(indx)+2*real(sqrt(varsigma{k}(indx)));flipdim(mu{k}(indx)-2*real(sqrt(varsigma{k}(indx))),1)];
    a = fill([Xtestsort; flipdim(Xtestsort,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(Xtestsort, Ytest{k}(indx), 'k--')];
    a =[ a plot(Xtestsort, mu{k}(indx),'k-')];
    c =plot(Xtrain{k},Ytrain{k},'ko', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
    minimum = min(mu{k}-2*real(sqrt(varsigma{k})));
    maximum = max(mu{k}+2*real(sqrt(varsigma{k})));
    if isfield(model, 'X_u') && ~isempty(model.X_u);
        b = plot(model.X_u, minimum -0.2*abs(minimum), 'kx');
        %        b = plot(model.X_u, 5, 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
    ylim = [0.5 5];
    ylabelText = ['Tide Height (m)'];
    ylabel(ylabelText,'fontsize',15)
    xlabelText = ['Time (days)'];
    xlabel(xlabelText,'fontsize',15)
    set(a,   'lineWidth', 2);
    set(c,   'markersize', 7);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'ylim', ylim, 'Color', 'none')
    box on
    if saveFigures==1
        if model.isSparse
            fileName = ['Tide_prediction_sp' approx num2str(k)];
        else
            fileName = ['Tide_prediction_full' num2str(k)];
        end
        print('-depsc', ['./results/' fileName]);
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
    %    zeroAxes(gca, [], 10, 'arial')
end

