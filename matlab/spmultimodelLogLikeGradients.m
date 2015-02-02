function [gParam, gX_u] = spmultimodelLogLikeGradients(model)

% SPMULTIMODELLOGLIKEGRADIENTS Compute the gradients for the parameters and X_u.
% FORMAT
% DESC
% Computes the gradients of the parameters in a spmultimodel structure. It
% uses sparsekernGradient.m to compute the gradients.
% ARG model : model for which the gradients are required.
% RETURN gParam : gradients of the parameters.
% RETURN gX_u : gradients of the inducing points.
%
% COPYRIGHT : Mauricio A Alvarez, 2009
%
% MODIFICATIONS : Mauricio A Alvarez, 2010

% MULTIGP

gX_u = [];


if model.varS
    initt = cputime;
    if model.isSpeed == 1 && isfield(model, 'subSpeed') && model.subSpeed == 2
        [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultimodelLocalCovGradient2(model);
        fhandle = str2func([model.betaTransform 'Transform']);
        gBeta = gBeta.*fhandle(model.beta, 'gradfact');
    else
        if isfield(model, 'beta') && ~isempty(model.beta)
            [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultimodelLocalCovGradient(model);
            fhandle = str2func([model.betaTransform 'Transform']);
            gBeta = gBeta.*fhandle(model.beta, 'gradfact');
        else
            [dLdKyy, dLdKuy, dLdKuu, dLdmu] = spmultimodelLocalCovGradient(model);
            gBeta = [];
        end
    end
    endtt = cputime - initt;
    if isfield(model, 'gamma') && ~isempty(model.gamma)
        gGamma = zeros(1, model.nlf);
        for i =1:model.nlf
            gGamma(i) = trace(dLdKuu{i});
        end
        fhandle = str2func([model.gammaTransform 'Transform']);
        gGamma = gGamma.*fhandle(model.gamma, 'gradfact');
    else
        gGamma = [];
    end
    if isfield(model, 'gammas') && ~isempty(model.gammas)
        gGammas = 0.5*(1./model.gammas - model.exS2(:)');
        fhandle = str2func([model.gammasTransform 'Transform']);
        gGammas = gGammas.*fhandle(model.gammas, 'gradfact');
    else
        gGammas = [];
    end
    switch model.isSpeed        
        case 1
            fhandle = str2func([model.kernType 'KernGradient']);
            gParam = fhandle(model.kern, model.outX, model.latX, dLdKyy, dLdKuy, dLdKuu);
        case {2,3}
            startVal = 1;
            endVal = 0;
            gParam = zeros(1, model.kern.nParams);
            submodel.nlf = 1;
            for i=1:model.nlf
                endVal = endVal + model.kern.comp{i}.nParams;
                tdLdKyy = dLdKyy(:,i);
                tdLdKuy = dLdKuy(i,:);
                tdLdKuu = dLdKuu(i);
                submodel.kern.comp{1} = model.kern.comp{i};
                submodel.kern.paramGroups = sparseDiag(model.kern.comp{i}.nParams);
                submodel.kernType = model.kernType;
                localX = cell(1+model.nout,1);
                localX{1,1} = model.latX{i};
                localX(2:end) = model.outX;
                submodel.X = localX;
                if model.isSpeed == 2
                    submodel.kernFuncNames = model.kernFuncNames;
                    gParam(startVal:endVal) = spmultimodelKernGradient(submodel, tdLdKyy, tdLdKuy, tdLdKuu);
                else
                    gParam(startVal:endVal) = sparseKernGradient(submodel, tdLdKyy, tdLdKuy, tdLdKuu);
                end
                startVal = endVal + 1;
            end
        otherwise
            error('Unknown SPEED option')
    end
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
    gParam = [gParam  g_meanFunc gGamma gBeta gGammas];
    if isfield(model, 'fix')
        for i = 1:length(model.fix)
            gParam(model.fix(i).index) = 0;
        end
    end
else
    switch model.approx
        case 'ftc'
            % Not implemented yet
        case {'fitc', 'pitc', 'dtcvar'}
            % Sparse approximations.
            if isfield(model, 'beta') && ~isempty(model.beta)
                [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultimodelLocalCovGradient(model);
                fhandle = str2func([model.betaTransform 'Transform']);
                gBeta = gBeta.*fhandle(model.beta, 'gradfact');
            else
                [dLdKyy, dLdKuy, dLdKuu, dLdmu] = spmultimodelLocalCovGradient(model);
                gBeta = [];
            end
            if isfield(model, 'gamma') && ~isempty(model.gamma)
                gGamma = zeros(1, model.nlf);
                for i =1:model.nlf
                    gGamma(i) = trace(dLdKuu{i});
                end
                fhandle = str2func([model.gammaTransform 'Transform']);
                gGamma = gGamma.*fhandle(model.gamma, 'gradfact');
            else
                gGamma = [];
            end
            startVal = 1;
            endVal = 0;
            gParam = zeros(1, model.kern.nParams);
            submodel.nlf = 1;
            for i=1:model.nlf
                endVal = endVal + model.kern.comp{i}.nParams;
                whichOutput = model.indOutGivenLat{i};
                tdLdKyy = dLdKyy(whichOutput);
                tdLdKuy = dLdKuy{i};
                tdLdKuu = dLdKuu(i);
                submodel.kern.comp{1} = model.kern.comp{i};
                submodel.kern.paramGroups = sparseDiag(model.kern.comp{i}.nParams);
                submodel.kernType = model.kernType;
                localX = cell(1+model.numOutGivenLat(i),1);
                localX{1,1} = model.latX{i};
                for j = 1:model.numOutGivenLat(i),
                    localX{1+j,1} = model.outX{whichOutput(j)};
                end
                submodel.X = localX;
                gParam(startVal:endVal) = sparseKernGradient(submodel, tdLdKyy, tdLdKuy, tdLdKuu);
                startVal = endVal + 1;
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
            gParam = [gParam  g_meanFunc gGamma gBeta];
            if isfield(model, 'fix')
                for i = 1:length(model.fix)
                    gParam(model.fix(i).index) = 0;
                end
            end
    end
end
