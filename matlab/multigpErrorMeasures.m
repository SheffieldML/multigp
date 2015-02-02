function [mae, mse, smse, msll, correl] = multigpErrorMeasures(y, yTest, mu, varsigma, nout, type)

% MULTIGPERRORMEASURES Measures SMSE and MSLL for multigp.
% FORMAT
% DESC computes SMSE and MSLL for the multi output Gaussian process. 
% ARG y : training data. A cell where each component is a training vector.
% The dimension is nout.
% ARG yTest : testing data. A cell where each component is a testing
% vector. The dimension is nout.
% ARG mu : predictive mean. The dimension is nout.
% ARG varsigma : predictive variance. The dimension is nout.
% ARG nout : number of outputs
% ARG type : if 'ind' considers each output as independent, if 'dep'
% considers dependency among otputs
% RETURN mae  : mean absolute error
% RETURN mse  : mean square error
% RETURN smse : standardized mean square error
% RETURN msll : mean standardized log loss
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if nargin < 6
    type = 'ind';
end

switch type
    case 'ind'
        mae = zeros(1, nout);
        mse = zeros(1, nout);
        smse = zeros(1, nout);
        msll = zeros(1, nout);
        snlp = zeros(1, nout);
        correl = zeros(1, nout);
        for k = 1:nout,
            maerror = abs(yTest{k} - mu{k});
            mserror = (yTest{k} - mu{k}).^2;
            std_mserror = mserror/var(yTest{k});
            sserror = ((yTest{k} - mu{k}).^2)./(2*varsigma{k});
            logprob = 0.5*log(2*pi.*varsigma{k}) + sserror - 0.5*log(2*pi.*var(y{k}))...
                - ((yTest{k} - mean(y{k})).^2)./(2*var(y{k}));
            mae(k) = mean(maerror);
            mse(k) = mean(mserror);
            smse(k) = mean(std_mserror);
            msll(k) = mean(logprob);
            snlp(k)= errorMeasureRegress(mu{k}, varsigma{k}, yTest{k}, 'snlp', y{k});
            cR2 = corrcoef(yTest{k}, mu{k});
            correl(k) = cR2(2,1)^2;
        end
    case 'dep'
        yVec = cell2mat(y);
        yTestVec = cell2mat(yTest);
        muVec = cell2mat(mu);
        varsigmaVec = cell2mat(varsigma);
        globalGaussian = true;        
        if globalGaussian
            meanyTrain = mean(yVec);
            varyTrain = var(yVec);           
        else
            meanyTrain = zeros(length(yTestVec),1);
            varyTrain = zeros(length(yTestVec),1);
            startOne = 1;
            endOne = 0;
            for i=1:nout
                endOne = endOne + length(yTest{i});
                meanyTrain(startOne:endOne) = mean(y{i});
                varyTrain(startOne:endOne) = var(y{i});
                startOne = endOne + 1;
            end
        end
        maerror = abs(yTestVec - muVec);
        mserror = (yTestVec - muVec).^2;
        std_mserror = mserror/var(yTestVec);
        sserror = ((yTestVec - muVec).^2)./(2*varsigmaVec);
        logprob = 0.5*log(2*pi.*varsigmaVec) + sserror - 0.5*log(2*pi*varyTrain)...
            - ((yTestVec - meanyTrain).^2)./(2*varyTrain);
        mae = mean(maerror);
        mse = mean(mserror);
        smse = mean(std_mserror);
        msll = mean(logprob);
        snlp = errorMeasureRegress(muVec, varsigmaVec, yTestVec, 'snlp', yVec);
end


% With Michalis version
% for k = 1:nout,
%     mse(k) = errorMeasureRegress(mu{k}, varsigma{k}, yTest{k}, 'mse', y{k});
%     smse(k) = errorMeasureRegress(mu{k}, varsigma{k}, yTest{k}, 'smse', y{k});
%     msll(k) = errorMeasureRegress(mu{k}, varsigma{k}, yTest{k}, 'snlp', y{k});
% end






