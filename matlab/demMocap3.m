% DEMMOCAP3 Demo of Sparse Multi Output Gaussian Process using PITC and the lfm kernel 
% FORMAT
% DESC 


% MULTIGP

% clc
clear
rand('seed',1e6);
randn('seed',1e6);

% load('./results/paramInitFull7Out');
% rand('seed',seed);
% randn('seed',seed);

%
% dataSetName = 'datasetODE31_30_5Ind1';
% dataSetName = 'datasetMOCAPSubject49Test016Dec32SelChan7';
% dataSetName = 'datasetMOCAPSubject49Test016Dec32SelChan2';
% dataSetName = 'datasetMOCAPSubject49Train18Test19Dec8-32SelChan9-Balance-RightArm';
dataSetName = 'datasetMOCAPSubject49Train18Test19Dec32SelChan56-Balance-Full';
experimentNo = 1;
ntrainX =200;
ntrainX2 =15;
iters =3000;
approx = 'pitc';
saveFigures=1;
% load data
data = multigpLoadData(dataSetName);
data.nin = 1;  % Number of latent functions
missingData = cell(data.nout,1);

options.isSparse = 1;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 1;  % Indicates if the outputs are to be standarized
options.outKernelName = 'sim'; % Indicates the name of the output kernel: lfm, sim or gg
options.tieSelectMethod = 'free';


% Setup model
[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsMocap(options, data, ntrainX, missingData, ntrainX2, approx);

options.optimiser = 'conjgrad';
options.includeInd = 1;
% Creates the model
model = multigpCreate(MYtrain, MXtrain, options);


% Change initial parameters

lengthScaleConv = 0.5*(max(data.Xtrain{1}) + rand(model.nlf,1));
%lengthScaleConv = 0.5*Delta_t;
inverseWidthConv = 1./(lengthScaleConv.^2);
%inverseWidthConv = 1./(lengthScaleConv);
for k =1: model.nlf
    paramsReal.(['b' num2str(k)]) = inverseWidthConv(k);
end
if options.includeInd % length scale of the size of the interpoint distance
    lengthScaleInd = (data.Xtrain{1}(2) -  data.Xtrain{1}(1)) + rand(model.nout,1);
%     lengthScaleInd = delta_t*ones(model.nout,1);
%    lengthScaleInd = rand(model.nout,1);
    inverseWidthInd = 1./(lengthScaleInd.^2);
    for k =1:model.nout
        paramsReal.(['l' num2str(k)]) = inverseWidthInd(k);
    end
end
paramsChanged = multigpChangeParamInit( options, paramsReal, ...
    model);
model= multigpExpandParam(model, paramsChanged);

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
        b = plot(model.X_u, minimum-abs(0.5*minimum), 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
%    ylim = [minimum-abs(0.5*minimum) maximum+0.2*maximum];
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

%Zt = linspace(min(Xtrain{1})-0.2,max(Xtrain{1})+0.2,100)';
Zt = linspace(0,12,200)';
% This part is to plot the posterior distribution
[mean_u, varsigma_u] = multigpPosteriorMeanVar(model, Zt);
% xlim = [min(model.X_u) max(model.X_u)];
for j =1:model.nlf
    figure    
    hold on
    f = [mean_u{j}+2*real(sqrt(varsigma_u{j}));flipdim(mean_u{j}-2*real(sqrt(varsigma_u{j})),1)];
    a = fill([Zt; flipdim(Zt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
%     a =[ a plot(data.Ztrain{j}, data.Utrain{j}, 'k--')];
    a =[ a plot(Zt,mean_u{j},'k-')];
    set(a,   'lineWidth', 2);
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none')
    if saveFigures==1
        fileName = ['Toy_posterior' approx num2str(j) '_outputs_' num2str(options.nout)];
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
