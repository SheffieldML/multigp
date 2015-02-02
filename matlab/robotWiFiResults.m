function meanAbsError = robotWiFiResults(model, XTest, yTest)

% MULTIGP




if strcmp(model.type, 'multimodel')
    meanAbsError = zeros(model.numModels,1);
    for k=1:model.numModels,
        muP = gpPosteriorMeanVar(model.comp{k}, XTest{k});
        meanAbsError(k) = mean(abs(muP - yTest{k}));
    end
else
    meanAbsError = zeros(model.nout,1);
    for k=1:model.nout,
        if strcmp(model.approx, 'ftc')
            XtestMod  = XTest{k};
        else
            XtestMod = XTest{k};
        end
        Ytest = yTest{k};
        muP = multigpPosteriorMeanVar(model,  XtestMod);
        kk = muP{model.nlf+k};
        meanAbsError(k) = mean(abs(Ytest - kk));
    end
end