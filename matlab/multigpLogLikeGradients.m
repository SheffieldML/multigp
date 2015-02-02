function [gParam, gX_u] = multigpLogLikeGradients(model)

% MULTIGPLOGLIKEGRADIENTS Compute the gradients for the parameters and X_u.

% COPYRIGHT : Mauricio A Alvarez, 2008, 2010

% MULTIGP

gX_u = [];
g_scaleBias = gpScaleBiasGradient(model);

if isfield(model, 'isSpeedUp') && ~isempty(model.isSpeedUp) && (model.isSpeedUp == 2)
    if strcmp(model.approx, 'ftc')
        covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
        covGrad = 0.5*covGrad;
        gParam = globalKernGradient(model.kern, model.X(model.nlf+1:end), covGrad);
        if ~model.includeNoise && isfield(model, 'beta') && ~isempty(model.beta)
            gBeta = zeros(1, model.nout);
            startVal = 1;
            endVal = 0;
            for k=1:model.nout
                endVal = endVal + size(model.X{model.nlf+k},1);
                gBeta(k) = - model.beta(k)^(-2)*trace(covGrad(startVal:endVal,startVal:endVal));
                startVal = endVal + 1;
            end
            fhandle = str2func([model.betaTransform 'Transform']);
            gBeta = gBeta.*fhandle(model.beta, 'gradfact');
        else
            gBeta = [];
        end
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                    || strcmp(model.kernType,'lfmwhite') ...
                    || strcmp(model.kernType,'simwhite') ...
                    || strcmp(model.kernType,'simnorm') ...
                    || strcmp(model.kernType,'simglobal')
                gmuFull = model.m'*model.invK;
                gmu = zeros(1, model.nout);
                startVal = 1;
                endVal = 0;
                for j=1:model.nout,
                    endVal =  endVal + size(model.X{model.nlf+j},1);
                    gmu(j) = sum(gmuFull(startVal:endVal));
                    startVal = endVal + 1;
                end
                g_meanFunc = meanGradient(model.meanFunction, gmu);
            else
                g_meanFunc = [];
            end
        else
            g_meanFunc = [];
        end
        gParam = [gParam g_scaleBias g_meanFunc gBeta];
    else
        % Sparse approximations. It does not provide yet gradients with respect
        % to the inducing variables.
        if isfield(model, 'beta') && ~isempty(model.beta)
            [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultigpLocalCovGradient(model);
            fhandle = str2func([model.betaTransform 'Transform']);
            gBeta = gBeta.*fhandle(model.beta, 'gradfact');
        else
            [dLdKyy, dLdKuy, dLdKuu, dLdmu] = spmultigpLocalCovGradient(model);
            gBeta = [];
        end
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            gGamma = zeros(1, model.nlf);
            if strcmp(model.kernType, 'lmc')
                for i =1:model.nlf
                    for j=1:model.rankCorregMatrix
                        linI = (i-1)*model.rankCorregMatrix + j;
                        gGamma(i) = gGamma(i) + trace(dLdKuu{linI});
                    end
                end
            else
                for i =1:model.nlf
                    gGamma(i) = trace(dLdKuu{i});
                end
            end
            fhandle = str2func([model.gammaTransform 'Transform']);
            gGamma = gGamma.*fhandle(model.gamma, 'gradfact');
        else
            gGamma = [];
        end
        gParam = globalKernGradient(model.kern, model.X(model.nlf+1:end), ...
            dLdKyy,  model.X(1:model.nlf), dLdKuy, dLdKuu);
%         if ~model.fixInducing
%             [gParam, gX_u] = sparseKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
%         else
%             gParam = sparseKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
%         end
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                    || strcmp(model.kernType,'lfmwhite') ...
                    || strcmp(model.kernType,'simwhite') ...
                    || strcmp(model.kernType,'simnorm') ...
                    || strcmp(model.kernType,'simglobal')
                gmu = zeros(1,model.nout);
                for j=1:model.nout,
                    gmu(j) = sum(dLdmu{j});
                end
                g_meanFunc = meanGradient(model.meanFunction, gmu);
            else
                g_meanFunc = [];
            end
        else
            g_meanFunc = [];
        end
        gParam = [gParam g_scaleBias g_meanFunc gGamma gBeta];
        if isfield(model, 'fix')
            for i = 1:length(model.fix)
                gParam(model.fix(i).index) = 0;
            end
        end
        % if there is only one output argument, pack gX_u and gParam into it.
        if ~model.fixInducing || nargout > 1
            gParam = [gParam gX_u(:)'];
        end
    end
else
    switch model.approx
        case 'ftc'
            covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
            covGrad = 0.5*covGrad;
            if strcmp(model.kernType, 'lmc')
                gParam = zeros(1, model.kern.nParams);
                startVal = 1;
                endVal = 0;
                for i=1:length(model.kern.comp)
                    for j=1:model.kern.comp{i}.numBlocks
                        endVal = endVal + model.kern.comp{i}.comp{j}.nParams;
                        if j>model.nlf
                            if model.isIsotopic
                                gParam(startVal:endVal) = kernGradient(model.kern.comp{i}.comp{j}, ...
                                    model.X{j}, covGrad);
                            else
                                gParam(startVal:endVal) = kernGradient(model.kern.comp{i}.comp{j}, ...
                                    model.X(j:end), covGrad);
                            end
                        end
                        startVal = endVal + 1;
                    end
                end
                base = 0;
            else
                if sum(cellfun('length',model.X)) == model.N,
                    base = 0;
                else
                    base = model.nlf;
                end
                covGrad2 = zeros(base + model.N);
                covGrad2(base+1:base+model.N,base+1:base+model.N) = covGrad;
                covGrad = covGrad2;
                gParam = kernGradient(model.kern, model.X, covGrad);
            end
            if ~model.includeNoise && isfield(model, 'beta') && ~isempty(model.beta)
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && model.noiseOpt
                    covGrad = covGrad(base+1:base+model.N,base+1:base+model.N);
                    gBeta = zeros(1, model.nout);
                    startVal = 1;
                    endVal = 0;
                    for k=1:model.nout
                        endVal = endVal + size(model.X{model.nlf+k},1);
                        gBeta(k) = - model.beta(k)^(-2)*trace(diag((1./model.nRepeats{k}))...
                            *covGrad(startVal:endVal,startVal:endVal));
                        startVal = endVal + 1;
                    end
                    fhandle = str2func([model.betaTransform 'Transform']);
                    gBeta = gBeta.*fhandle(model.beta, 'gradfact');
                end
            else
                gBeta = [];
            end
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                        || strcmp(model.kernType,'lfmwhite') ...
                        || strcmp(model.kernType,'simwhite') ...
                        || strcmp(model.kernType,'simnorm')
                    gmuFull = model.m'*model.invK;
                    gmu = zeros(1,model.d);
                    startVal = 1;
                    endVal = 0;
                    for j=1:model.d-base,
                        endVal =  endVal + length(model.X{j+base});
                        gmu(j+base) = sum(gmuFull(startVal:endVal));
                        startVal = endVal + 1;
                    end
                    gmu = gmu(model.nlf+1:end);
                    g_meanFunc = meanGradient(model.meanFunction, gmu);
                else
                    g_meanFunc = [];
                end
            else
                g_meanFunc = [];
            end
            gParam = [gParam g_scaleBias g_meanFunc gBeta];
            if isfield(model, 'fix')
                for i = 1:length(model.fix)
                    gParam(model.fix(i).index) = 0;
                end
            end
        case {'dtc', 'fitc', 'pitc', 'dtcvar'}
            % Sparse approximations.
            if isfield(model, 'beta') && ~isempty(model.beta)
                [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultigpLocalCovGradient(model);
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    if model.noiseOpt == 3
                        gBeta = [];
                    else
                        fhandle = str2func([model.betaTransform 'Transform']);
                        gBeta = gBeta.*fhandle(model.beta, 'gradfact');
                    end
                else
                    fhandle = str2func([model.betaTransform 'Transform']);
                    gBeta = gBeta.*fhandle(model.beta, 'gradfact');
                end
            else
                [dLdKyy, dLdKuy, dLdKuu, dLdmu] = spmultigpLocalCovGradient(model);
                gBeta = [];
            end
            if isfield(model, 'gamma') && ~isempty(model.gamma)
                gGamma = zeros(1, model.nlf);
                if strcmp(model.kernType, 'lmc')
                    for i =1:model.nlf
                        for j=1:model.rankCorregMatrix
                            linI = (i-1)*model.rankCorregMatrix + j;
                            gGamma(i) = gGamma(i) + trace(dLdKuu{linI});
                        end
                    end
                else
                    for i =1:model.nlf
                        gGamma(i) = trace(dLdKuu{i});
                    end
                end
                fhandle = str2func([model.gammaTransform 'Transform']);
                gGamma = gGamma.*fhandle(model.gamma, 'gradfact');
            else
                gGamma = [];
            end
            if ~model.fixInducing
                [gParam, gX_u] = sparseKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
            else
                gParam = sparseKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
            end
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                        || strcmp(model.kernType,'lfmwhite') ...
                        || strcmp(model.kernType,'simwhite') ...
                        || strcmp(model.kernType,'simnorm')
                    gmu = zeros(1,model.nout);
                    for j=1:model.nout,
                        gmu(j) = sum(dLdmu{j});
                    end
                    g_meanFunc = meanGradient(model.meanFunction, gmu);
                else
                    g_meanFunc = [];
                end
            else
                g_meanFunc = [];
            end
            gParam = [gParam g_scaleBias g_meanFunc gGamma gBeta];
            if isfield(model, 'fix')
                for i = 1:length(model.fix)
                    gParam(model.fix(i).index) = 0;
                end
            end
            % if there is only one output argument, pack gX_u and gParam into it.
            if ~model.fixInducing || nargout > 1
                gParam = [gParam gX_u(:)'];
            end
    end
end



