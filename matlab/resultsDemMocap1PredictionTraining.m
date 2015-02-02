function model = resultsDemMocap1PredictionTraining(dirResults, ...
    experimentName, subject, numberChannels, numberMissing, ...
    numberLatentFunctions, experiment, useInd);

% RESULTSDEMMOCAP1PREDICTIONTRAINING description

% MULTIGP

  %
% Plot Mocap database results

% Loading the data

if (experiment == 0)
    if useInd
        experimentName = strcat(dirResults, 'Mocap_Subject', num2str(subject), ...
            '_', experimentName, '_', num2str(numberChannels), 'outputs_', ...
            num2str(numberLatentFunctions), 'lf_Ind');
    else
        experimentName = strcat(dirResults, 'Mocap_Subject', num2str(subject), ...
            '_', experimentName, '_', num2str(numberChannels), 'outputs_', ...
            num2str(numberLatentFunctions), 'lf_NoInd');
    end
else
    if useInd
        experimentName = strcat(dirResults, 'Mocap_Subject', num2str(subject), ...
            '_', experimentName, '_', num2str(numberChannels), 'outputs_', ...
            num2str(numberLatentFunctions), 'lf_Ind_Exp', num2str(experiment));
    else
        experimentName = strcat(dirResults, 'Mocap_Subject', num2str(subject), ...
            '_', experimentName, '_', num2str(numberChannels), 'outputs_', ...
            num2str(numberLatentFunctions), 'lf_NoInd_Exp', num2str(experiment));
    end
end
eval(['load ' experimentName])

% Code required to use models trained old version of multigpCreate

if (model.type == 'lfm')
    model.type = 'multigp';
    model.kernType = 'lfm';
end

% Display initial length-scale

if exist('invWidthlfm', 'var')
    disp(' ');
    disp(['Initial length-scale = ' num2str(1./sqrt(invWidthlfm))]);
end

error_aprox = zeros(numberChannels, 1);

disp(' ');
disp('MSE    SNR(dB)');
disp(' ');

saveFigures = 1;

model.Ytest = Ytest;
model.Ytrain = Ytrain;
model.Xtest = Xtest;
for i=1:model.nout
    options.biasYteststd(i) = std(Ytest{i});
    options.biasYtestmean(i) = mean(Ytest{i});
end

% This part is to show the plots of the mean prediction and error bars

Xt = linspace(min(Xtrain{1})-0.2,max(Xtrain{1})+0.2,300)';
% Xt = Xtrain{1};
[mu, varsigma] = multigpPredictionMeanVarTraining(model, Xt, numberMissing, options);
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
    yr = spline(data.Xtrain{k}, data.Ytrain{k}, Xt);
    energy_yr(k) = mean(yr.*yr);
    error_approx(k) = mean((mu{k}-yr).^2);
    disp(['Channel ' num2str(k) ': ' num2str([error_approx(k) 10*log10(energy_yr(k)/error_approx(k))])])
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
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', [Zt(1) Zt(end)], 'Color', 'none')
    box on
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

disp(' ');
disp(['Average  : ' num2str([mean(error_approx) 10*log10(mean(energy_yr./error_approx))])]);

disp(' ');
disp(['Minus Log Likelihood = ' num2str(modelLogLikelihood(model))]);

disp(' ');
disp(['Time (min.) = ' num2str(total_t/60)]);

return
