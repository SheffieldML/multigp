function CC = schoolResults(model, XTest, yTest)

% SCHOOLRESULTS description. 
  
% MULTIGP
  
CC = zeros(1, model.nout);

for k= 1:model.nout,
    if strcmp(model.approx, 'ftc')
        XtestMod  = XTest{k};
    else
        XtestMod = XTest{k};
    end
    Ytest = yTest{k};
    muP = multigpPosteriorMeanVar(model,  XtestMod);
    kk = muP{model.nlf+k};
    cR2 = corrcoef(Ytest,kk);
    CC(k) = cR2(2,1)^2; 
end
