function [Y, fPos] = multigpOut(model, x, f, varargin)

% MULTIGPOUT outputs of a multigp given the latent forces
% FORMAT
% DESC evaluates the outputs of a multigp Gaussian process model for a set
% of inputs x
% ARG model : the model for which the output is being evaluated.
% ARG x : the input position for which the output is required.
% RETURN y : the output of the MULTIGP model.
%
% FORMAT
% DESC evaluates the outputs of a multigp Gaussian process model for a set
% of inputs x and values of latent functions given
% ARG model : the model for which the output is being evaluated.
% ARG x : the input position for which the output is required.
% ARG f : latent function value for which the output is required.
% RETURN y : the output of the MULTIGP model.
%
% SEEALSO : multigpCreate, multigpPosteriorMeanVar
%
% COPYRIGHT : Mauricio Alvarez, 2009

% MULTIGP

sampling = false;

if isfield(model, 'scaleVal') && ~isempty(model.scaleVal)
    model.scaleVal = model.scaleVal/4;
else
    model.scaleVal = 1;
end

% if sampling
%     rand('twister', 1e6);
%     randn('state', 1e6);
% end
%

if ~isempty(varargin)
   initPos = varargin{1};
end

if nargin < 3,
    % Here it should just do multigpPosteriorMeanVar
else
    if ~iscell(x)
        X = cell(model.d,1);
        for i = 1:model.d
            X{i} = x;
        end
    end
    if length(f) ~= model.nlf,
        error('The number of latent functions in f must be equal to the number of latent functions in the model');
    end
    switch model.approx,
        case {'ftc'}
            Kyf = kernCompute(model.kern, x);
            dim = length(x) ;
            Kffinvf = zeros(dim*model.nlf,1);
            startVal = 1;
            endVal = 0;
            for k = 1: model.nlf,
                endVal = endVal + dim;
                [Kffinv, U, jitter] = pdinv(Kyf(startVal:endVal,startVal:endVal));
                if jitter > 1e-4
                    fprintf('Warning: multigpOutput added jitter of %2.4f\n', jitter)
                end
                Kffinvf(startVal:endVal,1) = Kffinv*f{k};
                startVal = endVal + 1;
            end
            fullMean = Kyf(startVal:end,1:dim*model.nlf)*Kffinvf;            
            if sampling
                covar =  Kyf(startVal:end,startVal:end);                
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, X, model.nlf);
                end
                m = m(model.nlf*dim + 1:end,1);
                startVal=1;
                endVal=0;                
                for i=1:model.nout
                    endVal = endVal + dim;
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        fullMean(startVal:endVal, 1)  = (fullMean(startVal:endVal, 1)*model.scale(i+model.nlf) ...
                            + m(startVal:endVal, 1)+ model.bias(i+model.nlf))*model.scaleVal;
                    else
                        fullMean(startVal:endVal, 1)  = (fullMean(startVal:endVal, 1)*model.scale(i+ model.nlf) + model.bias(i+model.nlf))*...
                            model.scaleVal;
                    end
                    startVal = endVal+1;
                end                
                y = real(gsamp(fullMean, covar*model.scaleVal^2, 1));                
                y = y';
                startVal=1;
                endVal=0;                               
                Y = zeros(dim, model.nout);
                for j=1:model.nout,
                    endVal = endVal + dim;
                    Y(:,j) = y(startVal:endVal, 1);
                    startVal = endVal+1;
                end
            else
                y = fullMean;
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, X, model.nlf);
                end
                startVal=1;
                endVal=0;
                m = m(model.nlf*dim + 1:end,1);
                Y = zeros(dim, model.nout);
                for i=1:model.nout
                    endVal = endVal + dim;
                    if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                        Y(:,i)  = (y(startVal:endVal, 1)*model.scale(i+model.nlf) + m(startVal:endVal, 1)...
                            + model.bias(i+model.nlf))*model.scaleVal;
                    else
                        Y(:,i)  = (y(startVal:endVal, 1)*model.scale(i+model.nlf) + model.bias(i+model.nlf))*...
                            model.scaleVal;
                    end
                    startVal = endVal+1;
                end
            end

        case {'dtc', 'fitc', 'pitc','dtcvar'}
            dim = length(x);
            if strcmp(model.kernType, 'lfmglobal')
                kernTemp = model.kern;
                kernTemp.approx = 'dtc';
                [~, Kyu, Kuu] = globalKernCompute(kernTemp, X, X(1:model.nlf));
                Kuuinv =  cell(model.nlf,1);
                Kuuinvu =  cell(model.nlf,1);
                for r = 1:model.nlf
                    [Kuuinv{r}, ~, jitter] = pdinv(Kuu{r});
                    if ~sampling
                        Kuuinvu{r} = Kuuinv{r}*f{r};
                    end
                    
                end
            else                
                Kuu =  cell(model.nlf,1);
                Kuuinv =  cell(model.nlf,1);
                sqrtKuu = cell(model.nlf,1);
                Kuuinvu = cell(model.nlf,1);
                Kyu =  cell(model.nout, model.nlf);
                Kyy =  cell(model.nout,1);
                y = cell(model.nout,1);
                % Compute Kuu
                for r = 1:model.nlf,
                    Kuu{r,1} = zeros(dim);
                    for c = 1: length(model.kern.comp)
                        Kuu{r,1} = Kuu{r,1} + real(multiKernComputeBlock(model.kern.comp{c},  x, r, r));
                    end
                    [Kuuinv{r}, sqrtKuu{r}, jitter] = pdinv(Kuu{r});
                    if ~sampling
                        Kuuinvu{r} = Kuuinv{r}*f{r};
                    end
                    if jitter > 1e-4
                        fprintf('Warning: multigpOutput added jitter of %2.4f\n', jitter)
                    end
                    for i =1:model.nout,
                        Kyu{i,r} = real(multiKernComputeBlock(model.kern.comp{r},  x,...
                            x, i+model.nlf, r));
                    end
                end
            end
            % Compute the Kyy part
            if sampling
                %ys = cell(model.nout,1);
                Y = zeros(dim, model.nout);
                KyuKuuinv = cell(model.nout, model.nlf);
                KyuKuuinvKuy = cell(model.nout, 1);                
                D = cell(model.nout, 1);
                Dinv = cell(model.nout, 1);
                KuyDinv = cell(model.nlf, model.nout);
                KuyDinvy = cell(model.nlf, 1);
                for i=1:size(KuyDinvy, 1),
                    KuyDinvy{i} = zeros(dim,1);
                end
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
                end
                for j =1:model.nout,
                    switch model.approx
                        case {'dtc','fitc'}
                            Kyy{j,1} = zeros(dim,1);
                        case 'pitc'
                            Kyy{j,1} = zeros(dim);
                    end
                    for c = 1: length(model.kern.comp)
                        switch model.approx
                            case {'dtc','fitc'}
                                Kyy{j,1} = Kyy{j,1} + real(kernDiagCompute(model.kern.comp{c}.comp{j+model.nlf}, x));
                            case 'pitc'
                                Kyy{j,1} = Kyy{j,1} + real(multiKernComputeBlock(model.kern.comp{c}, x, j+model.nlf, j+model.nlf));
                        end
                    end
                    for r =1:model.nlf,
                        KyuKuuinv{j,r} = Kyu{j,r}*Kuuinv{r};
                    end
                    KyuKuuinvKuy{j} = zeros(dim);
                    y{j} = zeros(dim,1);
                    for k =1:model.nlf,
                        KyuKuuinvKuy{j} = KyuKuuinvKuy{j} + KyuKuuinv{j,k}*Kyu{j,k}';
                        y{j,1} = y{j,1} + KyuKuuinv{j,r}*f{r};
                    end
                    % endVal = endVal + dim;
                    switch model.kernType
                        case {'gg','ggwhite'}
                            y{j,1} = y{j,1}*model.scale(j) + model.bias(j);
                        otherwise
                            y{j,1} = y{j,1}*model.scale(j) + m(j) + model.bias(j);
                    end
                    y{j,1} = y{j,1}*model.scaleVal;
                    if nargout>1
                        D{j} = Kyy{j} - KyuKuuinvKuy{j};
                        Dinv{j} = pdinv(D{j});
                    end
                    Y(:,j) = real(gsamp(y{j}, (Kyy{j} - KyuKuuinvKuy{j})*(model.scaleVal^2), 1));
                end
                
                if nargout > 1

                    for k =1:model.nlf,
                        for q =1:model.nout,
                            switch model.approx
                                case 'dtc'
                                    %
                                case {'fitc','pitc'}
                                    KuyDinv{k,q} = Kyu{q,k}'*Dinv{q};
                            end
                        end
                    end

                    for r =1:model.nlf,
                        KuyDinvy{r,1} = zeros(dim,1);
                        for q =1:model.nout,
                            KuyDinvy{r} = KuyDinvy{r} + KuyDinv{r,q}*Y(:,q);
                        end
                    end

                    A = cell(model.nlf);

                    for k =1:model.nlf,
                        for r =1:model.nlf,
                            KuyDinvKyu = zeros(dim);
                            for q =1:model.nout,
                                KuyDinvKyu = KuyDinvKyu + KuyDinv{k,q}*Kyu{q,r};
                            end
                            if (k == r)
                                A{k,r} = Kuu{k} + KuyDinvKyu;
                            else
                                A{k,r} = KuyDinvKyu;
                            end
                        end
                    end

                    AinvMat = pdinv(cell2mat(A));
                    Ainv  = mat2cell(AinvMat,dim*ones(1,model.nlf), dim*ones(1,model.nlf));
                    AinvKuyDinvy = cell(model.nlf,1);

                    for r = 1:model.nlf,
                        AinvKuyDinvy{r,1} = zeros(dim,1);
                        for k = 1:model.nlf,
                            AinvKuyDinvy{r} = AinvKuyDinvy{r} + Ainv{r,k}*KuyDinvy{k};
                        end
                    end

                    % Mean of the posterior
                    fPos = cell(model.nlf, 1);
                    for r = 1:model.nlf
                        fPos{r} = Kuu{r}*AinvKuyDinvy{r}/model.scaleVal;
                    end

                end

            else
                Y = zeros(dim, model.nout);
%                 if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
%                     m1 = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
%                 end
                m = initPos/model.scaleVal; 
                outputs = zeros(model.nout*dim,1);
                startVal = 1;
                endVal = 0;
                for q =1: model.nout,
                    endVal = endVal + dim;
                    for k =1:model.nlf,
                        outputs(startVal:endVal) = outputs(startVal:endVal) + Kyu{q,r}*Kuuinvu{r};
                    end
                    switch model.kernType
                        case {'gg','ggwhite'}
                            Y(:,q) = outputs(startVal:endVal)*model.scale(q) + model.bias(q);
                        otherwise
%                             Y(:,q) = outputs(startVal:endVal)*model.scale(q) + model.bias(q);
                             Y(:,q) = outputs(startVal:endVal)*model.scale(q) + m(q) + model.bias(q);
                    end
                    Y(:,q) = Y(:,q)*model.scaleVal;
                    startVal = endVal + 1;
                end
            end

        otherwise
    end
end
