function [mu, varsig] = spmultimodelPosteriorMeanVar(model, X, computeAll)

% SPMULTIMODELPOSTERIORMEANVAR gives mean and variance of the posterior distribution.
% FORMAT
% DESC gives the mean and variance of outputs for the sparse multimodel 
% ARG model : the model for which posterior is to be computed.
% ARG X : cell array containing locations where outputs are to be computed.
% RETURN mu : cell array containing mean posterior vectors.
% RETURN varsig : cell array containing the variance posterior vector
%
% COPYRIGHT : Mauricio A Alvarez, 2008, 2009, 2010

% MULTIGP

% If computeAll is true, then compute posterior for latent forces and
% outputs. If computeAll is false, then compute posterior only for latent
% forces. This is introduced particularly for applications in which the
% number of outputs >> number of latent forces and we are only interested
% in latent forces posteriors.
if nargin<3
    computeAll = true;
end

% If X is a vector assume it applies for all outputs.
if ~iscell(X)
    xtemp = X;
    X = cell(1, model.d);
    for i = 1:model.d
        X{i} = xtemp;
    end
end

if model.varS
    % Computing the full posterior in this case is intractable, so we only
    % compute the mean function and the covariance function associated to
    % each output
    if computeAll
        mu = cell(model.nout+model.nlf,1);
        varsig = cell(model.nout+model.nlf,1);
    else
        mu = cell(model.nlf,1);
        varsig = cell(model.nlf,1);
    end
    Ku_star_u    = cell(model.nlf,1);
    KX_star_X2   = cell(model.nout,1);
    KX_star_Diag = cell(model.nout,model.nlf);
    % Precomputations
    switch model.isSpeed
        case 1
            fhandle = str2func([model.kernType 'KernComputeTest']);
            if isfield(model, 'gamma') && ~isempty(model.gamma)
                [KX_star_Diag, KX_star_X2, Ku_star_u] = fhandle(model.kern, model.latX, X, X{1}, model.gamma);
            else
                [KX_star_Diag, KX_star_X2, Ku_star_u] = fhandle(model.kern, model.latX, X, X{1});
            end
        case {2,3}            
            for r = 1:model.nlf,
                Ku_star_u{r,1} = real(multiKernComputeBlock(model.kern.comp{r}, X{1}, model.latX{r}, 1, 1));
                if isfield(model, 'gamma') && ~isempty(model.gamma)
                    Ku_star_u{r,1} =  Ku_star_u{r,1} + model.gamma(r)*eye(size(Ku_star_u{r,1}));
                end
            end
            for i = 1:model.nout,
                for j = 1: model.nlf
                    % Compute Kff
                    KX_star_Diag{i,j} = real(kernDiagCompute(model.kern.comp{j}.comp{1+i}, X{i}));
                    KX_star_X2{i,j} = real(multiKernComputeBlock(model.kern.comp{j},  X{i},...
                        model.latX{j}, i+1, 1));
                end
            end
    end
    % Posterior for the latent functions
    for j=1:model.nlf
        mu{j} = Ku_star_u{j}*model.AinvKuyHatDinvy{j};
        varsig{j} = diag(Ku_star_u{j}*model.Ainv{j,j}*Ku_star_u{j}');% This is because we are only interested in the variances
    end
    % Posterior for the output functions
    if computeAll       
        KuuinvAinv = cell(model.nlf);
        for r =1:model.nlf,
            for q =1:model.nlf,
                if r ==q,
                    KuuinvAinv{r,q} = model.Kuuinv{r} - model.Ainv{r,r} ...
                        + model.AinvKuyHatDinvy{r}*model.AinvKuyHatDinvy{r}';
                else
                    KuuinvAinv{r,q} = -model.Ainv{r,q} ...
                        + model.AinvKuyHatDinvy{r}*model.AinvKuyHatDinvy{q}';
                end
            end
        end
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
        end
        for i=1:model.nout,
            mux_star = zeros(size(X{i},1),1);            
            Kx_starXu = zeros(size(X{i},1),1);
            Kx_star_x_star = zeros(size(X{i},1),1); 
            %%% Faster Alternative
            for r =1:model.nlf,
                mux_star = mux_star + model.qs.mean(i,r)*...
                    KX_star_X2{i,r}*model.AinvKuyHatDinvy{r};
                exS2 = model.qs.Sigma(r,r,i) + model.qs.mean(i,r)^2;
                temp = KuuinvAinv{r,r}*KX_star_X2{i,r}';
                Kx_star_x_star = Kx_star_x_star + exS2*(KX_star_Diag{i,r} ...
                    - sum(KX_star_X2{i,r}.*temp',2));
                for q =1:r-1,
                    exS2 = model.qs.Sigma(r,q,i) + model.qs.mean(i,r)*model.qs.mean(i,q);
                    temp = KuuinvAinv{r,q}*KX_star_X2{i,q}';
                    Kx_starXu = Kx_starXu + 2*exS2*sum(KX_star_X2{i,r}.*temp',2);
                end
            end
            if model.includeInd
                % Not implemented yet.
            end
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                mu{i+model.nlf} = mux_star*model.scale(i) + m(i) + model.bias(i);
            else
                mu{i+model.nlf} = mux_star*model.scale(i) + model.bias(i);
            end            
            if nargout == 2
                varsig{i+model.nlf} = (Kx_star_x_star -  Kx_starXu ...
                    + 1/model.beta(i) +  mux_star.*mux_star)*model.scale(i)*model.scale(i);
            end
        end
    end
else
    switch model.approx
        case 'ftc'
            % Not implemented yet
        case {'fitc','pitc','dtcvar'}
            if computeAll
                mu = cell(model.nout+model.nlf,1);
                varsig = cell(model.nout+model.nlf,1);
            else
                mu = cell(model.nlf,1);
                varsig = cell(model.nlf,1);
            end
            Ku_star_u    = cell(model.nlf,1);
            KX_star_X2   = cell(model.nout,1);
            KX_star_Diag = cell(model.nout,1);
            latGivenOut = zeros(1, model.nout);
            % Precomputations
            for r = 1:model.nlf,
                Ku_star_u{r,1} = real(multiKernComputeBlock(model.kern.comp{r}, X{1}, model.latX{r}, 1, 1));
                if isfield(model, 'gamma') && ~isempty(model.gamma)
                    Ku_star_u{r,1} =  Ku_star_u{r,1} + model.gamma(r)*eye(size(Ku_star_u{r,1}));
                end
                whichOutputs = model.indOutGivenLat{r};
                latGivenOut(whichOutputs) = latGivenOut(whichOutputs) + 1;
                for j =1:model.numOutGivenLat(r);
                    temp = real(multiKernComputeBlock(model.kern.comp{r},  X{whichOutputs(j)},...
                        model.latX{r}, j+1, 1));
                    KX_star_X2{whichOutputs(j),1}{latGivenOut(whichOutputs(j)),1} = temp;
                end
            end
            % Posterior for the latent functions
            for j=1:model.nlf
                mu{j} = Ku_star_u{j}*model.AinvKuyDinvy{j};
                varsig{j} = diag(Ku_star_u{j}*model.Ainv{j,j}*Ku_star_u{j}');% This is because we are only interested in the variances
            end
            % Posterior for the output functions
            if computeAll
                % Compute the Kyy part
                for j =1:model.nout,
                    latentFunc = model.indLatGivenOut{j};
                    KX_star_Diag{j} = zeros(size(X{j},1),1);
                    for i = 1: length(latentFunc)
                        whichOutputs = model.indOutGivenLat{latentFunc(i)};
                        localOutput = find(j == whichOutputs );
                        KX_star_Diag{j} = KX_star_Diag{j} +  ...
                            real(kernDiagCompute(model.kern.comp{latentFunc(i)}.comp{1+localOutput}, X{j}));
                    end
                end
                KuuinvAinv = cell(model.nlf);
                for r =1:model.nlf,
                    for q =1:model.nlf,
                        if r ==q,
                            KuuinvAinv{r,q} = model.Kuuinv{r} - model.Ainv{r,r};
                        else
                            KuuinvAinv{r,q} = -model.Ainv{r,q};
                        end
                    end
                end
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
                end
                for i=1:model.nout,
                    DinvKyuAinv = cell(model.numLatGivenOut(i),1);
                    DinvKyuAinvKuyDinvy = zeros(model.sizeX(i),1);
                    mux_star = zeros(size(X{i},1),1);
                    Kx_starXu = zeros(size(X{i},1));
                    whichLatent = model.indLatGivenOut{i};
                    for r =1:model.numLatGivenOut(i),
                        whichOutputs1 = model.indOutGivenLat{whichLatent(r)};
                        subsetLatent1 = find(whichOutputs1 == i);
                        DinvKyuAinv{r} = zeros(model.sizeX(i),model.k);
                        for q =1:model.numLatGivenOut(i),
                            whichOutputs2 = model.indOutGivenLat{whichLatent(q)};
                            subsetLatent2 = find(whichOutputs2 == i);
                            Kx_starXu = Kx_starXu + KX_star_X2{i}{r}*...
                                KuuinvAinv{whichLatent(r),whichLatent(q)}*KX_star_X2{i}{q}';
                            DinvKyuAinv{r} = DinvKyuAinv{r} + model.KuyDinv{whichLatent(q)}{subsetLatent2}'...
                                *model.Ainv{whichLatent(q),whichLatent(r)};
                        end
                        mux_star = mux_star + KX_star_X2{i}{r}*...
                            model.AinvKuyDinvy{whichLatent(r)};
                        DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + model.KuyDinv{whichLatent(r)}{subsetLatent1}'*...
                            model.AinvKuyDinvy{whichLatent(r)};
                    end
                    DinvKyuAinvKuyDinv = zeros(model.sizeX(i));
                    DinvKyuAinvKuyStar = zeros(model.sizeX(i), size(X{i},1));
                    for r =1:model.numLatGivenOut(i),
                        whichOutputs1 = model.indOutGivenLat{whichLatent(r)};
                        subsetLatent1 = find(whichOutputs1 == i);
                        DinvKyuAinvKuyDinv = DinvKyuAinvKuyDinv +  DinvKyuAinv{r}*model.KuyDinv{whichLatent(r)}{subsetLatent1};
                        DinvKyuAinvKuyStar = DinvKyuAinvKuyStar +  DinvKyuAinv{r}*KX_star_X2{i}{r}';
                    end
                    if model.includeInd
                        % Not implemented yet.
                    end
                    switch model.kernType
                        case {'gg','ggwhite'}
                            mu{i+model.nlf} = mux_star*model.scale(i) + model.bias(i);
                        otherwise
                            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                                mu{i+model.nlf} = mux_star*model.scale(i) + m(i) + model.bias(i);
                            else
                                mu{i+model.nlf} = mux_star*model.scale(i) + model.bias(i);
                            end
                    end
                    if nargout == 2
                        if model.includeInd
                            if isfield(model, 'beta') && ~isempty(model.beta)
                                varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) - diag(covInd) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                            else
                                varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) - diag(covInd))*model.scale(i)*model.scale(i);
                            end
                        else
                            if isfield(model, 'beta') && ~isempty(model.beta)
                                varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                            else
                                varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu))*model.scale(i)*model.scale(i);
                            end
                        end
                    end
                end
            end
        otherwise
            error('spmultimodelPosteriorMeanVar not yet implemented');
    end
end