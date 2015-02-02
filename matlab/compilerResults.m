function [totalError, elapsed_time] = compilerResults(model, XTest, yTest)

% COMPILERRESULTS description.
  
% MULTIGP

maxTest = length(yTest{1});
numPart = 1000;
step = floor(maxTest/numPart);
if strcmp(model.approx, 'ftc')
    XtestMod = cell(model.nout+ model.nlf,1);
    for j=1:model.nlf
        XtestMod{j} = ones(1,13);
    end
else
    XtestMod = cell(model.nout,1);
end
mserror = zeros(maxTest,1);

elapsed_time = 0;
for j =1:(numPart+1),
    if (j==numPart+1)
        indexes = indexes(end)+1:maxTest;
    else
        indexes = (j-1)*step+1:j*step;
    end
    for k= 1:model.nout,
        if strcmp(model.approx, 'ftc')
            XtestMod{k+model.nlf} = XTest{k}(indexes,:); 
        else
            XtestMod{k} = XTest{k}(indexes,:);
        end
    end
    tic
    [mu, void] = multigpPosteriorMeanVar(model, XtestMod);
    localTime = toc;
    elapsed_time = elapsed_time + localTime; 
    for k= 1:model.nout,
        mserror(indexes,k) = abs((yTest{k}(indexes) - mu{k+model.nlf}));
    end
end

totalError = mean(mserror)';