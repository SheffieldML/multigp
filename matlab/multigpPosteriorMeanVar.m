function [mu, varsig, covar] = multigpPosteriorMeanVar(model, X, computeAll)

% MULTIGPPOSTERIORMEANVAR gives mean and variance of the posterior distribution.
% FORMAT
% DESC gives the mean and variance of outputs for the multi-output
% Gaussian process.
% ARG model : the model for which posterior is to be computed.
% ARG X : cell array containing locations where outputs are to be computed.
% RETURN mu : cell array containing mean posterior vectors.
% RETURN varsig : cell array containing the variance posterior vector
% RETURN covar : cell array containing the covariance posterior only for
% the latent forces
%
% DESC gives the mean and variance of outputs for the multi-output
% Gaussian process.
% ARG model : the model for which posterior is to be computed.
% ARG X : cell array containing locations where outputs are to be computed.
% ARG computeAll : a flag to indicate if mean posterior for outputs are to
% be computed. It is true by default.
% RETURN mu : cell array containing mean posterior vectors.
% RETURN varsig : cell array containing the variance posterior vector
% RETURN covar : cell array containing the covariance posterior only for
% the latent forces
%
% COPYRIGHT : Mauricio A Alvarez, 2008, 2009, 2010
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010
%
% SEEALSO : gpPosteriorMeanVar

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

if isfield(model, 'isSpeedUp') && ~isempty(model.isSpeedUp) && (model.isSpeedUp == 2)
    switch model.approx
        case 'ftc'
            mu = cell(model.d,1);
            varsig = cell(model.d,1);
            KX_star = globalKernCompute(model.kern, {model.X(model.nlf+1:end), X(model.nlf+1:end)});
            muTemp = KX_star'*model.alpha;
            diagK = globalKernDiagCompute(model.kern, X(model.nlf+1:end));
            diagKPlusNoise = diagK;
            if ~model.includeNoise && isfield(model, 'beta') && ~isempty(model.beta)
                startVal = 1;
                endVal = 0;
                for k=1:model.nout
                    endVal = endVal + size(X{model.nlf+k},1);
                    diagKPlusNoise(startVal:endVal,1) = diagKPlusNoise(startVal:endVal,1) + (1/model.beta(k));
                    startVal = endVal + 1;
                end
            end
            Kinvk = model.invK*KX_star;
            varsigTemp = diagKPlusNoise - sum(KX_star.*Kinvk, 1)';
            % Compute the mean and variance of the latent force
            % Hack the kernel and use it to compute the required
            % covariances
            kernTemp = model.kern;
            kernTemp.approx = 'dtc';
            [Kyy, Kyu, Kuu] = globalKernCompute(kernTemp, model.X(model.nlf+1:end), X(model.nlf+1:end));
            muForce = cell(model.nlf,1);            
            if nargout >2 && ~computeAll
                varForce = cell(model.nlf);
                for i =1:model.nlf
                    KuyLocal = cell2mat(Kyu(:,i))';
                    muForce{i} = KuyLocal*model.alpha;
                    varForce{i,i} = diag(Kuu{i}) - sum(KuyLocal.*temp',2);
                    for j =1:i-1
                        KuyLocal2 = cell2mat(Kyu(:,j))';
                        temp = model.invK*KuyLocal2';
                        varForce{i,j} = - sum(KuyLocal.*temp',2);
                        varForce{j,i} = (varForce{i,j})';
                    end
                end
                for i =1:model.nlf
                    mu{i} = muForce{i};
                    varsig{i} = varForce(i,:);
                end                
            else                
                varForce = cell(model.nlf,1);
                for i =1:model.nlf
                    KuyLocal = cell2mat(Kyu(:,i))';
                    muForce{i} = KuyLocal*model.alpha;
                    temp = model.invK*KuyLocal';
                    varForce{i} = diag(Kuu{i}) - sum(KuyLocal.*temp',2);
                end
            end
            if computeAll
                muTemp = [cell2mat(muForce); muTemp];
                varsigTemp = [cell2mat(varForce); varsigTemp];
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, X, model.nlf);
                end
                startVal=1;
                endVal=0;
                for i=1:length(X)
                    endVal = endVal + size(X{i},1);
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        mu{i}  = muTemp(startVal:endVal, 1)*model.scale(i) ...
                            + m(startVal:endVal, 1)+ model.bias(i);
                    else
                        mu{i}  = muTemp(startVal:endVal, 1)*model.scale(i) + model.bias(i);
                    end
                    varsig{i} = varsigTemp(startVal:endVal, 1)*model.scale(i)*model.scale(i);
                    startVal = endVal+1;
                end
            end
        case {'dtc', 'fitc', 'pitc', 'dtcvar'}
            if computeAll
                mu = cell(model.nout+model.nlf,1);
                varsig = cell(model.nout+model.nlf,1);
            else
                mu = cell(model.nlf,1);
                varsig = cell(model.nlf,1);
            end            
            kernTemp = model.kern;
            kernTemp.approx = 'dtc';
            [void, KX_star_X2, Ku_star_u] = globalKernCompute(kernTemp, ...
                X, {X(1:model.nlf), model.X(1:model.nlf)});
            KX_star_Diag = globalKernDiagCompute(model.kern, X);            
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
            % Posterior for the latent functions
            if nargout >2 && ~computeAll
                varForce = cell(model.nlf);
                for j=1:model.nlf
                    mu{j} = Ku_star_u{j}*model.AinvKuyDinvy{j};
                    varForce{j,j} = diag(Ku_star_u{j}*model.Ainv{j,j}*Ku_star_u{j}');% This is because we are only interested in the variances
                    for k=1:j-1
                        varForce{j,k} = diag(Ku_star_u{j}*model.Ainv{j,k}*Ku_star_u{k}');% This is because we are only interested in the variances
                        varForce{k,j} = (varForce{j,k})';
                    end
                end
                for j=1:model.nlf
                   varsig{j} = varForce(j,:);                   
                end
            else
                for j=1:model.nlf
                    mu{j} = Ku_star_u{j}*model.AinvKuyDinvy{j};
                    varsig{j} = diag(Ku_star_u{j}*model.Ainv{j,j}*Ku_star_u{j}');% This is because we are only interested in the variances
                end
            end
            % Posterior for the output functions
            if computeAll
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
                end                
                for i=1:model.nout,                    
                    mux_star = zeros(size(X{i},1),1);
                    Kx_starXu = zeros(size(X{i},1));
                    for r =1:model.nlf,                      
                        for q =1:model.nlf,
                            Kx_starXu = Kx_starXu + KX_star_X2{i,r}*KuuinvAinv{r,q}*KX_star_X2{i,q}';                         
                        end
                        mux_star = mux_star + KX_star_X2{i,r}*model.AinvKuyDinvy{r};                      
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
                        if isfield(model, 'beta') && ~isempty(model.beta)
                            varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                        else
                            varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu))*model.scale(i)*model.scale(i);
                        end
                    end
                end
            end
    end
else
    switch model.approx
        case 'ftc'
            mu = cell(model.d,1);
            varsig = cell(model.d,1);
            if strcmp(model.kernType, 'lmc')
                A = cell(1, model.nlf);
                B = cell(1, model.nlf); % Coregionalization matrices
                for r =1:model.nlf
                    A{r} = model.kern.comp{r}.comp{model.nlf+1}.A;
                    B{r} = model.kern.comp{r}.comp{model.nlf+1}.B;
                end
                KuuTemp = cell(1, model.nlf);
                KfuTemp = cell(1, model.nlf);
                for i=1:model.nlf
                    KuuTemp{i} = kernCompute(model.kern.comp{i}.comp{i}, X{i});
                    for l = 1:model.rankCorregMatrix
                        startOne = 1;
                        endOne = 0;
                        for j=1:model.nout
                            if model.isIsotopic
                                endOne = endOne + size(model.X{model.nlf+1},1);
                                KfuTemp{i,l}(startOne:endOne, :) = A{i}(j,l)*kernCompute(model.kern.comp{i}.comp{i}, ...
                                    model.X{model.nlf+1}, X{i});
                            else
                                endOne = endOne + size(model.X{model.nlf+j},1);
                                KfuTemp{i,l}(startOne:endOne, :) = A{i}(j,l)*kernCompute(model.kern.comp{i}.comp{i}, ...
                                    model.X{model.nlf+j}, X{i});
                            end
                            startOne = endOne + 1;
                        end
                    end
                end
                for i=1:model.nlf
                    for j=1:model.rankCorregMatrix
                        mu{i}{j} = KfuTemp{i,j}'*model.alpha;
                        varsig{i}{j} = diag(KuuTemp{i} - KfuTemp{i,j}'*model.invK*KfuTemp{i,j});
                    end
                end
                cont = model.nlf;
                KX_star = zeros(model.N, sum(cellfun('size', X(model.nlf+1:end),1)));
                diagK = zeros(sum(cellfun('size', X(model.nlf+1:end),1)),1);
                for i=1:length(model.kern.comp)
                    if model.isIsotopic
                        KX_star = KX_star + kernCompute(model.kern.comp{i}.comp{model.nlf+1}, ....
                            model.X{model.nlf+1}, X{model.nlf+1});
                    else
                        KX_star = KX_star + kernCompute(model.kern.comp{i}.comp{model.nlf+1}, ....
                            model.X(model.nlf+1:end), X(model.nlf+1:end));
                    end
                    diagK = diagK + kernDiagCompute(model.kern.comp{i}.comp{model.nlf+1}, ...
                        X(model.nlf+1:end));
                end
                muTemp = KX_star'*model.alpha;
                Kinvk = model.invK*KX_star;
                varsigTemp = diagK - sum(KX_star.*Kinvk, 1)';
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, X, model.nlf);
                end
                startVal=1;
                endVal=0;
                for i=1:model.nout
                    cont = cont + 1;
                    endVal = endVal + size(X{model.nlf+i},1);
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        mu{cont}  = muTemp(startVal:endVal, 1)*model.scale(model.nlf+i) ...
                            + m(startVal:endVal, 1)+ model.bias(model.nlf+i);
                    else
                        mu{cont}  = muTemp(startVal:endVal, 1)*model.scale(model.nlf+i) + model.bias(model.nlf+i);
                    end
                    varsig{cont} = varsigTemp(startVal:endVal, 1)*model.scale(model.nlf+i)*model.scale(model.nlf+i);
                    startVal = endVal+1;
                end
            else
                KX_star = kernCompute(model.kern, model.X, X);
                % Again, this is a bit of a hack to allow the inclusion of the
                % latent structure in the kernel structure.
                if sum(cellfun('length',model.X)) == model.N,
                    base = 0;
                else
                    base = model.nlf;
                end
                KX_star = KX_star(base+1:end,:);
                muTemp = KX_star'*model.alpha;
                diagK = kernDiagCompute(model.kern, X);
                Kinvk = model.invK*KX_star;
                varsigTemp = diagK - sum(KX_star.*Kinvk, 1)';
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, X, model.nlf);
                end
                startVal=1;
                endVal=0;
                for i=1:length(X)
                    endVal = endVal + size(X{i},1);
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        mu{i}  = muTemp(startVal:endVal, 1)*model.scale(i) ...
                            + m(startVal:endVal, 1)+ model.bias(i);
                    else
                        mu{i}  = muTemp(startVal:endVal, 1)*model.scale(i) + model.bias(i);
                    end
                    varsig{i} = varsigTemp(startVal:endVal, 1)*model.scale(i)*model.scale(i);
                    startVal = endVal+1;
                end
            end
        case {'dtc','fitc','pitc','dtcvar'}
            if strcmp(model.kernType, 'lmc')
                if computeAll
                    mu = cell(model.nout+model.nlf,1);
                    varsig = cell(model.nout+model.nlf,1);
                else
                    mu = cell(model.nlf,1);
                    varsig = cell(model.nlf,1);
                end
                A = cell(1, model.nlf);
                B = cell(1, model.nlf); % Coregionalization matrices
                for r =1:model.nlf
                    A{r} = model.kern.comp{r}.comp{model.nlf+1}.A;
                    B{r} = model.kern.comp{r}.comp{model.nlf+1}.B;
                end
                Ku_star_u    = cell(model.nlf,1);
                KX_star_X2   = cell(model.nout,model.nlf);
                for r = 1:model.nlf,
                    Ku_star_u{r,1} = multiKernComputeBlock(model.kern.comp{r}, X{1},  model.X{r}, r, r);
                    for i=1:model.nout
                        KX_star_X2{i,r} = real(multiKernComputeBlock(model.kern.comp{r}, X{i}, model.X{r}, r, r));
                        
                    end
                end
                KX_star_Diag = cell(model.nout,1);
                for i=1:model.nout
                    KX_star_Diag{i} = zeros(size(X{i},1),1);
                    for r=1:model.nlf
                        KX_star_Diag{i} = KX_star_Diag{i} + ...
                            B{r}(i,i)*kernDiagCompute(model.kern.comp{r}.comp{r}, X{i});
                    end
                end
                startOne = 1;
                endOne = 0;
                Ac = cell(model.nlf);
                KuuinvAinv = cell(model.nlf);
                for r=1:model.nlf
                    endOne = endOne + model.k(r)*model.rankCorregMatrix;
                    subAinvM = model.Ainv(startOne:endOne, startOne:endOne);
                    Ac{r,r} = mat2cell(subAinvM, model.k(r)*ones(model.rankCorregMatrix,1), ...
                        model.k(r)*ones(model.rankCorregMatrix,1));
                    for i=1:model.rankCorregMatrix
                        KuuinvAinv{r,r}{i,i} = model.Kuuinv{r} - Ac{r,r}{i,i};
                        for j=1:i-1
                            KuuinvAinv{r,r}{i,j} = - Ac{r,r}{i,j};
                            KuuinvAinv{r,r}{j,i} = - Ac{r,r}{j,i};
                        end
                    end
                    startTwo = 1;
                    endTwo = 0;
                    for k =1:r-1
                        endTwo = endTwo + model.k(k)*model.rankCorregMatrix;
                        subAinvM = model.Ainv(startOne:endOne, startTwo:endTwo);
                        Ac{r, k} = mat2cell(subAinvM, model.k(r)*ones(model.rankCorregMatrix,1), ...
                            model.k(k)*ones(model.rankCorregMatrix,1));
                        subAinvM = model.Ainv(startTwo:endTwo, startOne:endOne);
                        Ac{k, r} = mat2cell(subAinvM, model.k(k)*ones(model.rankCorregMatrix,1), ...
                            model.k(r)*ones(model.rankCorregMatrix,1));
                        startTwo = endTwo + 1;
                        for i=1:model.rankCorregMatrix
                            for j=1:model.rankCorregMatrix
                                KuuinvAinv{r,k}{i,j} = - Ac{r,k}{i,j};
                                KuuinvAinv{r,k}{j,i} = - Ac{r,k}{j,i};
                                KuuinvAinv{k,r}{i,j} = - Ac{k,r}{i,j};
                                KuuinvAinv{k,r}{j,i} = - Ac{k,r}{j,i};
                            end
                        end
                    end
                    startOne = endOne + 1;
                end
                % Posterior for the latent functions
                for j=1:model.nlf
                    for i=1:model.rankCorregMatrix
                        mu{j}{i} = Ku_star_u{j}*model.AinvKuyDinvy{j}(:,i);
                        varsig{j}{i} = diag(Ku_star_u{j}*Ac{j,j}{i,i}*Ku_star_u{j}');% This is because we are only interested in the variances
                    end
                end
                % Posterior for the output functions
                if computeAll
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
                    end
                    DinvKyuAinv = cell(model.nlf,1);
                    for i=1:model.nout,
                        DinvKyuAinvKuyDinvy = zeros(size(model.X{i+model.nlf},1),1);
                        mux_star = zeros(size(X{i},1),1);
                        Kx_starXu = zeros(size(X{i},1));
                        for r =1:model.nlf,
                            for n=1:model.rankCorregMatrix
                                DinvKyuAinv{r}{n} = zeros(size(model.X{i+model.nlf},1),model.k(r));
                                for q =1:model.nlf,
                                    for j=1:model.rankCorregMatrix
                                        Kx_starXu = Kx_starXu + ...
                                            A{r}(i,n)*KX_star_X2{i,r}*KuuinvAinv{r,q}{n,j}*(A{q}(i,j)*KX_star_X2{i,q}');
                                        DinvKyuAinv{r}{n} = DinvKyuAinv{r}{n} + ...
                                            (model.KuyDinv{q,i}(:,:,j))'*Ac{q,r}{j,n};
                                    end
                                end
                                mux_star = mux_star + A{r}(i,n)*KX_star_X2{i,r}*model.AinvKuyDinvy{r}(:,n);
                                DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                                    model.KuyDinv{r,i}(:,:,n)'*model.AinvKuyDinvy{r}(:,n);
                            end
                        end
                        DinvKyuAinvKuyDinv = zeros(size(model.X{i+model.nlf},1));
                        DinvKyuAinvKuyStar = zeros(size(model.X{i+model.nlf},1), size(X{i},1));
                        for r =1:model.nlf,
                            for j=1:model.rankCorregMatrix
                                DinvKyuAinvKuyDinv = DinvKyuAinvKuyDinv + ...
                                    DinvKyuAinv{r}{j}*(model.KuyDinv{r,i}(:,:,j));
                                DinvKyuAinvKuyStar = DinvKyuAinvKuyStar + ...
                                    DinvKyuAinv{r}{j}*(A{r}(i,j)*KX_star_X2{i,r}');
                            end
                        end
                        if model.includeInd
                            c = length(model.kern.comp);
                            if model.includeNoise
                                c = c - 1;
                            end
                            Kw_star_x2 = real(multiKernComputeBlock(model.kern.comp{c}, X{i}, model.X{i+model.nlf}, i+model.nlf, i+model.nlf));
                            Kw_star_fDinvKyuAinvKuyStar = Kw_star_x2*DinvKyuAinvKuyStar;
                            Kw_star_fKyyinvKfw_star = Kw_star_x2*(model.Dinv{i} - DinvKyuAinvKuyDinv)*Kw_star_x2';
                            covInd = Kw_star_fDinvKyuAinvKuyStar + Kw_star_fDinvKyuAinvKuyStar' + Kw_star_fKyyinvKfw_star;
                            muInd = Kw_star_x2*(model.Dinv{i}*model.m{i} - DinvKyuAinvKuyDinvy);
                            mux_star = mux_star + muInd;
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
                                    varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) ...
                                        - diag(covInd) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                                else
                                    varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) ...
                                        - diag(covInd))*model.scale(i)*model.scale(i);
                                end
                            else
                                if isfield(model, 'beta') && ~isempty(model.beta)
                                    varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) ...
                                        + 1/model.beta(i))*model.scale(i)*model.scale(i);
                                else
                                    varsig{i+model.nlf} = (KX_star_Diag{i} - ...
                                        diag(Kx_starXu))*model.scale(i)*model.scale(i);
                                end
                            end
                        end
                    end
                end
            else
                if computeAll
                    mu = cell(model.nout+model.nlf,1);
                    varsig = cell(model.nout+model.nlf,1);
                else
                    mu = cell(model.nlf,1);
                    varsig = cell(model.nlf,1);
                end
                Ku_star_u    = cell(model.nlf);
                KX_star_X2   = cell(model.nout,model.nlf);
                KX_star_Diag = cell(model.nout,1);
                for r = 1:model.nlf,
                    Ku_star_u{r,1} = zeros(size(X{1},1), model.k(r));
                    for c = 1: length(model.kern.comp)
                        Ku_star_u{r,1} = Ku_star_u{r,1} + ...
                            real(multiKernComputeBlock(model.kern.comp{c}, X{1}, model.X{r}, r, r));
                    end
                    for i =1:model.nout,
                        KX_star_X2{i,r} = real(multiKernComputeBlock(model.kern.comp{r},  X{i},...
                            model.X{r}, i+model.nlf, r));
                    end
                end
                % Compute the Kyy part
                for j =1:model.nout,
                    KX_star_Diag{j} = zeros(size(X{j},1),1);
                    for c = 1: length(model.kern.comp)
                        KX_star_Diag{j} = KX_star_Diag{j} + ...
                            real(kernDiagCompute(model.kern.comp{c}.comp{j+model.nlf}, X{j}));
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
                % Posterior for the latent functions
                for j=1:model.nlf
                    mu{j} = Ku_star_u{j}*model.AinvKuyDinvy{j};
                    varsig{j} = diag(Ku_star_u{j}*model.Ainv{j,j}*Ku_star_u{j}');% This is because we are only interested in the variances
                end
                % Posterior for the output functions
                if computeAll
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
                    end
                    DinvKyuAinv = cell(model.nlf,1);
                    for i=1:model.nout,
                        DinvKyuAinvKuyDinvy = zeros(size(model.X{i+model.nlf},1),1);
                        mux_star = zeros(size(X{i},1),1);
                        Kx_starXu = zeros(size(X{i},1));
                        for r =1:model.nlf,
                            DinvKyuAinv{r} = zeros(size(model.X{i+model.nlf},1),model.k(r));
                            for q =1:model.nlf,
                                Kx_starXu = Kx_starXu + KX_star_X2{i,r}*KuuinvAinv{r,q}*KX_star_X2{i,q}';
                                DinvKyuAinv{r} = DinvKyuAinv{r} + model.KuyDinv{q,i}'*model.Ainv{q,r};
                            end
                            mux_star = mux_star + KX_star_X2{i,r}*model.AinvKuyDinvy{r};
                            DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + model.KuyDinv{r,i}'*model.AinvKuyDinvy{r};
                        end
                        DinvKyuAinvKuyDinv = zeros(size(model.X{i+model.nlf},1));
                        DinvKyuAinvKuyStar = zeros(size(model.X{i+model.nlf},1), size(X{i},1));
                        for r =1:model.nlf,
                            DinvKyuAinvKuyDinv = DinvKyuAinvKuyDinv +  DinvKyuAinv{r}*model.KuyDinv{r,i};
                            DinvKyuAinvKuyStar = DinvKyuAinvKuyStar +  DinvKyuAinv{r}*KX_star_X2{i,r}';
                        end
                        if model.includeInd
                            c = length(model.kern.comp);
                            if model.includeNoise
                                c = c - 1;
                            end
                            Kw_star_x2 = real(multiKernComputeBlock(model.kern.comp{c}, X{i}, model.X{i+model.nlf}, i+model.nlf, i+model.nlf));
                            Kw_star_fDinvKyuAinvKuyStar = Kw_star_x2*DinvKyuAinvKuyStar;
                            Kw_star_fKyyinvKfw_star = Kw_star_x2*(model.Dinv{i} - DinvKyuAinvKuyDinv)*Kw_star_x2';
                            covInd = Kw_star_fDinvKyuAinvKuyStar + Kw_star_fDinvKyuAinvKuyStar' + Kw_star_fKyyinvKfw_star;
                            muInd = Kw_star_x2*(model.Dinv{i}*model.m{i} - DinvKyuAinvKuyDinvy);
                            mux_star = mux_star + muInd;
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
                                    if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                                        switch model.noiseOpt
                                            case 0
                                                varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                                            case 1
                                            case 2
                                            case 3
                                                varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) )*model.scale(i)*model.scale(i);
                                        end
                                    else
                                        varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                                    end
                                else
                                    varsig{i+model.nlf} = (KX_star_Diag{i} - diag(Kx_starXu))*model.scale(i)*model.scale(i);
                                end
                            end
                        end
                    end
                end
            end
        otherwise
            error('multigpPosteriorMeanVar not yet implemented');
    end
end



