% DEMGPSIM3_V2 Demo of full multi output GP with the sim kernel with PITC using data 2
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process with the sim kernel using
%      PITC and data 2

% MULTIGP

clc
clear
rand('seed',1e6);
randn('seed',1e6);
%
dataSetName = 'dataBarencoOption_0_Genes_5';
experimentNo = 1;
ntrainX =200;
ntrainX2 =5;
iters =3000;
approx = 'pitc';
saveFigures=1;
% load data
data = multigpLoadData(dataSetName);
for j=1:3,
    data.data{j}.nin = 2;  % Number of latent functions
end
%
options.isSparse = 1;      % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;      % Indicates if the inputs are to be standarized (1) or not(0)
options.stdBiasY = 0;      % Indicates if the outputs are to be standarized (1) or not(0)
options.outKernelName = 'sim'; % Indicates the name of the output kernel: lfm, sim or gg
options.inference = 'ml';  % Indicates the method of inference for sensitivities between ML and MAP
options.tieSelectMethod = 'within';

% Model type
model.type = 'cmultigp';
MXtrain = cell(3,1);
Xtrain = cell(3,1);
MYtrain = cell(3,1);
Ytrain = cell(3,1);
Ytest = cell(3,1);
Xtest = cell(3,1);

for j =1:3,
    % Setup each model
    [options, MXtrain{j}, Xtrain{j}, MYtrain{j}, Ytrain{j}, Xtest{j}, Ytest{j}] = multigpOptionsSim(options, data.data{j}, ntrainX,...
        ntrainX2, approx);
    options.optimiser = 'scg';
    % Creates the model
    model.comp{j} = multigpCreate(MYtrain{j}, MXtrain{j}, options);   
end


if options.nlf > 1,
    lengthScaleLf = [12 6];
    for k =1: model.comp{1}.nlf
        paramsReal.(['b' num2str(k)]) = 1/(lengthScaleLf(k)^2);
    end
    if options.includeInd % length scale of the size of the interpoint distance
        lengthScaleInd = ones(model.comp{1}.nout,1);
        inverseWidthInd = 1./(lengthScaleInd.^2);
        for k =1:model.comp{1}.nout
            paramsReal.(['l' num2str(k)]) = inverseWidthInd(k);
        end
    end
    for j =1:3
        paramsChanged = multigpChangeParamInit( options, paramsReal, ...
            model.comp{j} );
        model.comp{j} = multigpExpandParam(model.comp{j}, paramsChanged);
    end
end

model.optimiser =  'scg';

% Trains the model and counts the training time
ini_t = cputime;
model = modelOptimise(model, [], [], 1, iters);
total_t = cputime - ini_t;

%This part is to show the plots of the mean prediction and error bars
close all
if data.data{1}.nout<10,
    %for j =1:length(model.comp)
    for j =1:10
        Xt = linspace(min(Xtrain{j}{1})-0.2,max(Xtrain{j}{1})+0.2,100)';
        [mu, varsigma] = multigpPredictionMeanVar(model.comp{j}, Xt, options);
        xlim = [min(Xtrain{j}{1}) max(Xtrain{j}{1})];
        for k=1:model.comp{1}.nout,
            figure
            hold on
            f = [mu{k}+2*real(sqrt(varsigma{k}));flipdim(mu{k}-2*real(sqrt(varsigma{k})),1)];
            a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
            a =[ a plot(Xt, mu{k},'k-')];
            c =plot(Xtrain{j}{k},Ytrain{j}{k},'k.');
            minimum = min(mu{k}-2*real(sqrt(varsigma{k})));
            maximum = max(mu{k}+2*real(sqrt(varsigma{k})));
            if isfield(model, 'X_u') && ~isempty(model.X_u);
                b = plot(model.X_u, -8, 'kx');
                set(b, 'linewidth', 2)
                set(b, 'markersize', 10);
            end
            ylim = [min(minimum,min(Ytrain{j}{k})) max(maximum,max(Ytrain{j}{k}))];
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
    end
end

for i=1:length(model.comp)
    Zt = linspace(min(Xtrain{i}{1})-0.2,max(Xtrain{i}{1})+0.2,100)';
    % This part is to plot the posterior distribution
    [mean_u, varsigma_u] = multigpPosteriorMeanVar(model.comp{i}, Zt);
    % xlim = [min(model.X_u) max(model.X_u)];
    for j =1:model.comp{1}.nlf
        figure
        hold on
        f = [mean_u{j}+2*real(sqrt(varsigma_u{j}));flipdim(mean_u{j}-2*real(sqrt(varsigma_u{j})),1)];
        a = fill([Zt; flipdim(Zt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
        %   a =[ a plot(data.Ztrain{j}, data.Utrain{j}, 'k--')];
        a =[ a plot(Zt,mean_u{j},'k-')];
        set(a,   'lineWidth', 2);
        set(gca, 'fontname', 'arial', 'fontsize', 15, 'Color', 'none')
        if saveFigures==1
            fileName = ['Toy_posterior_' approx '_nlf_' num2str(model.comp{1}.nlf)   '_latentNumber_' num2str(j) ...
                '_Replica_' num2str(i) '_Z_' num2str(ntrainX2)];
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
end

save(['./results/modelSim_' approx '_NumLf' num2str(model.comp{1}.nout) 'NumNout' num2str(model.comp{1}.nout) ...
    'NumZ' num2str(ntrainX2)])

