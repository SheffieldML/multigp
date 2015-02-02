function model = spmultigpUpdateAD(model)

% SPMULTIGPUPDATEAD Update the representations of A and D associated with
% the model.
% FORMAT
% DESC updates the representations of A and D in the model when
% called by spmultigpUpdateKernels.
% ARG model : the model for which the representations are being
% updated.
% RETURN model : the model with the A and D representations
% updated.
%
% SEEALSO : spmultigpUpdateKernels, spmultigpExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010


% MULTIGP

logDetKuu = 0;
for r=1:model.nlf,
    % Kuuinv really is the inverse of Kuu + regularization term
    if isfield(model, 'gamma') && ~isempty(model.gamma)
        [model.Kuuinv{r}, model.sqrtKuu{r}] = pdinv(model.KuuGamma{r});
        model.logDetKuu{r} = logdet(model.KuuGamma{r}, model.sqrtKuu{r});
    else
        [model.Kuuinv{r}, model.sqrtKuu{r}] = pdinv(model.Kuu{r});
        model.logDetKuu{r} = logdet(model.Kuu{r}, model.sqrtKuu{r});
    end
    logDetKuu = logDetKuu + model.logDetKuu{r};
end
model.logDetKuuT = logDetKuu;
for r =1: model.nlf,
    for k =1: model.nout,
        model.KuuinvKuy{r,k} = model.Kuuinv{r}*model.Kyu{k,r}';
    end
end
logDetD = 0;
traceDinvyy = 0;
model.KtildeT = 0;


if strcmp(model.kernType, 'lmc')
    % Correction for the corregionalization model in the log of the
    % determinants of the inverse
    model.logDetKuuT = model.rankCorregMatrix*model.logDetKuuT;
    A = cell(1, model.nlf);
    B = cell(1, model.nlf); % Coregionalization matrices
    for r =1:model.nlf
        A{r} = model.kern.comp{r}.comp{model.nlf+1}.A;
        B{r} = model.kern.comp{r}.comp{model.nlf+1}.B;
    end
    switch model.approx
        case {'dtc', 'dtcvar'}
            for k=1:model.nout
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                        (model.noiseOpt == 1)
                    % This option is to use only in the School data
                    if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                        model.D{k} = sparseDiag( (1/model.beta(k))*(1./model.nRepeats{k}));
                        model.Dinv{k} = sparseDiag(model.beta(k)*model.nRepeats{k});
                        model.logDetD{k} = sum(log((1/model.beta(k))*(1./model.nRepeats{k})));
                        if strcmp(model.approx, 'dtcvar')
                            model.Ktilde{k} = zeros(size(model.Kyy{k},1),1);
                            for r =1: model.nlf,
                                model.Ktilde{k} = model.Ktilde{k} + B{r}(k,k)*(model.Kyy{k,r} - ...
                                    sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2));
                            end
                            temp = sum(model.beta(k)*model.nRepeats{k}.*model.Ktilde{k});
                        end
                    else
                        error('Model does not contain nRepeats')
                    end
                else
                    model.D{k} = sparseDiag(1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1));
                    model.Dinv{k} = sparseDiag(model.beta(k)*ones(size(model.X{k+model.nlf},1),1));
                    model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
                    if strcmp(model.approx, 'dtcvar')
                        model.Ktilde{k} = zeros(size(model.Kyy{k},1),1);
                        for r =1: model.nlf,
                            model.Ktilde{k} = model.Ktilde{k} + B{r}(k,k)*(model.Kyy{k,r} - ...
                                sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2));
                        end
                        temp = model.beta(k)*sum(model.Ktilde{k});
                    end
                end
                if strcmp(model.approx, 'dtcvar')
                    model.KtildeT = model.KtildeT + temp;
                end
                model.Dinvy{k} = model.Dinv{k}*model.m{k};
                logDetD = logDetD + model.logDetD{k};
                traceDinvyy = traceDinvyy + sum(diag(model.Dinv{k}).*model.m{k}.*model.m{k});
            end
            
        case 'fitc'
            for k=1:model.nout
                model.D{k}  = zeros(size(model.Kyy{k},1),1);
                for r =1: model.nlf,
                    model.D{k} = model.D{k} + B{r}(k,k)*(model.Kyy{k,r} - ...
                        sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2));
                end
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                        (model.noiseOpt == 1)
                    if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                        model.D{k} = (model.D{k} + (1/model.beta(k))*(1./model.nRepeats{k}));
                    else
                        error('Model does not contain nRepeats')
                    end
                else
                    model.D{k} = (model.D{k} + 1/model.beta(k));
                end
                %                 model.Dinv{k} = sparseDiag(1./model.D{k});
                %                 model.logDetD{k} = sum(log(model.D{k}));
                %                 model.D{k} = sparseDiag(model.D{k});% This is to keep the flow of the gradients
                %                 model.Dinvy{k} = model.Dinv{k}*model.m{k};
                %                 logDetD = logDetD + model.logDetD{k};
                %                 traceDinvyy = traceDinvyy + sum((model.Dinv{k}*model.m{k}).*model.m{k});
                %%%%%
                model.Dinv{k} = 1./model.D{k};
                model.logDetD{k} = sum(log(model.D{k}));
                model.Dinvy{k} = model.Dinv{k}.*model.m{k};
                logDetD = logDetD + model.logDetD{k};
                traceDinvyy = traceDinvyy + sum(model.Dinvy{k}.*model.m{k});
            end
        case 'pitc'
            for k=1:model.nout
                model.D{k}  = zeros(size(model.Kyy{k},1));
                for r =1: model.nlf,
                    model.D{k} = model.D{k} + B{r}(k,k)*(model.Kyy{k,r} - model.Kyu{k,r}*model.KuuinvKuy{r,k});
                end
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                        (model.noiseOpt == 1)
                    if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                        model.D{k} =(model.D{k} + diag((1/model.beta(k))*(1./model.nRepeats{k})));
                    else
                        error('Model does not contain nRepeats')
                    end
                else
                    model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
                end
                [model.Dinv{k}, model.sqrtD{k}] = pdinv(model.D{k});
                model.logDetD{k} = logdet(model.D{k}, model.sqrtD{k});
                model.Dinvy{k} = model.Dinv{k}*model.m{k};
                logDetD = logDetD + model.logDetD{k};
                %traceDinvyy = traceDinvyy + sum((model.Dinv{k}*model.m{k}).*model.m{k});
                traceDinvyy = traceDinvyy + sum(model.Dinvy{k}.*model.m{k});
            end
    end
    model.logDetDT = logDetD;
    model.traceDinvyy = traceDinvyy;
    
    for r =1:model.nlf,
        model.KuyDinvy{r,1} = zeros(model.k(r),model.rankCorregMatrix);
        for i=1:model.rankCorregMatrix
            for q =1:model.nout,
                model.KuyDinvy{r}(:,i) = model.KuyDinvy{r}(:,i) + ...
                    A{r}(q,i)*(model.Kyu{q,r}')*model.Dinvy{q};
                model.KuyDinv{r,q}(:,:,i) = A{r}(q, i)*(model.Kyu{q,r}')*model.Dinv{q};
            end
        end
    end
    
    Apos = cell(model.nlf*model.rankCorregMatrix);
    KuyDinvy = cell(model.nlf, 1);
    startValOne = 1;
    endValOne = 0;
    for k =1:model.nlf
        for i=1:model.rankCorregMatrix
            KuyDinvKyu = zeros(model.k(k));
            for q =1:model.nout,
                KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}(:,:,i)*(A{k}(q,i)*model.Kyu{q,k});
            end
            if isfield(model, 'gamma') && ~isempty(model.gamma)
                model.A{k,k}{i,i} = model.KuuGamma{k} + KuyDinvKyu;
            else
                model.A{k,k}{i,i} = model.Kuu{k} + KuyDinvKyu;
            end
            for j =1:i-1,
                KuyDinvKyu = zeros(model.k(k));
                for q =1:model.nout,
                    KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}(:,:,i)*(A{k}(q,j)*model.Kyu{q,k});
                end
                model.A{k,k}{i,j} = KuyDinvKyu;
                model.A{k,k}{j,i} = KuyDinvKyu';
            end
        end
        endValOne = endValOne + model.rankCorregMatrix;
        Apos(startValOne:endValOne, startValOne:endValOne) = model.A{k,k};
        startValTwo = 1;
        endValTwo = 0;
        for r =1:k-1,
            endValTwo = endValTwo + model.rankCorregMatrix;
            for i=1:model.rankCorregMatrix
                for j=1:model.rankCorregMatrix
                    KuyDinvKyu = zeros(model.k(k), model.k(r));
                    for q =1:model.nout,
                        KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}(:,:,i)*(A{r}(q,j)*model.Kyu{q,r});
                    end
                    model.A{k,r}{i,j} = KuyDinvKyu;
                    model.A{r,k}{j,i} = KuyDinvKyu';
                end
            end
            Apos(startValOne:endValOne, startValTwo:endValTwo) = model.A{k,r};
            Apos(startValTwo:endValTwo, startValOne:endValOne) = model.A{r,k};
            startValTwo = endValTwo + 1;
        end
        KuyDinvy{k} = model.KuyDinvy{k}(:);
        startValOne = endValOne + 1;
    end
    
    AM = cell2mat(Apos);
    [AMinv, sqrtAM, jitter]  = pdinv(AM);
    sqrtAMinv = sqrtAM\eye(size(sqrtAM, 1));
    KuyDinvyM = cell2mat(KuyDinvy);
    model.logDetA = logdet(AM, sqrtAM);
    model.Ainv = AMinv;
    model.sqrtAinvKuyDinvy = mat2cell((sqrtAMinv')*KuyDinvyM, cellfun('length', KuyDinvy), 1);
    AMinvKuyDinvy = AMinv*KuyDinvyM;
    
    startValOne = 1;
    endValOne = 0;
    for r = 1:model.nlf,
        for i=1:model.rankCorregMatrix
            endValOne = endValOne + model.k(r);
            model.AinvKuyDinvy{r}(:,i) = AMinvKuyDinvy(startValOne:endValOne);
            startValOne = endValOne + 1;
        end
    end
    model.AinvKuyDinvy2 = mat2cell(AMinvKuyDinvy, cellfun('length', KuyDinvy), 1);
else
    
    switch model.approx
        case 'dtc'
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0                        
                        for k =1:model.nout
                            model.D{k} = 1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
                            model.Dinv{k} = model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
                            model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
                            for r=1:model.nlf
                                model.KuyDinv{r,k} = model.beta(k)*model.Kyu{k,r}';
                            end
                            logDetD = logDetD + model.logDetD{k};
                            traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                        end
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            for k =1:model.nout
                                model.D{k} = (1/model.beta(k))*(1./model.nRepeats{k});
                                model.Dinv{k} = model.beta(k)*model.nRepeats{k};
                                model.logDetD{k} = sum(log((1/model.beta(k))*(1./model.nRepeats{k})));
                                for r=1:model.nlf
                                    tempo = model.Dinv{k}(:, ones(1, model.k(r)));
                                    model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';                                   
                                end
                                logDetD = logDetD + model.logDetD{k};
                                traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                            end
                        else
                            error('Model does not contain nRepeats')
                        end
                    case 2
                        for k =1:model.nout
                            % Be aware: the interpretation of beta for this option is as
                            % variance, not precision
                            model.D{k} = 1./model.beta{k};
                            model.Dinv{k} = model.beta{k};
                            model.logDetD{k} = log(1./model.beta{k});
                            for r=1:model.nlf
                                tempo = model.beta{k}(:, ones(1, model.k(r)));
                                model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';                   
                            end
                            logDetD = logDetD + model.logDetD{k};
                            traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                        end
                end
            else
                for k =1:model.nout
                    model.D{k}(1:size(model.X{k+model.nlf},1),1) = 1/model.beta(k);
                    model.Dinv{k}(1:size(model.X{k+model.nlf},1),1) = model.beta(k);
                    model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
                    for r=1:model.nlf
                        model.KuyDinv{r,k} = model.beta(k)*model.Kyu{k,r}';
                    end
                    logDetD = logDetD + model.logDetD{k};
                    traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                end
            end            
        case 'fitc'
            for k =1:model.nout
                KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
                for r =1: model.nlf,
                    KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2);
                end
                model.D{k} = model.Kyy{k} - KyuKuuinvKuy; % In fitc D is a diagonal matrix.
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            model.D{k} = (model.D{k} + 1/model.beta(k));
                        case 1
                            if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                                model.D{k} = (model.D{k} + (1/model.beta(k))*(1./model.nRepeats{k}));
                            else
                                error('Model does not contain nRepeats')
                            end
                        case 2
                            % Be aware: the interpretation of beta for this option is as
                            % variance, not precision
                            model.D{k} = (model.D{k} + 1./model.beta{k});
                    end
                else
                    model.D{k} = (model.D{k} + 1/model.beta(k));
                end
                %                 model.Dinv{k} = sparseDiag(1./model.D{k});
                %                 model.logDetD{k} = sum(log(model.D{k}));
                %                 model.D{k} = sparseDiag(model.D{k});% This is to keep the flow of the gradients
                %                 for r=1:model.nlf
                %                     model.KuyDinv{r,k} = model.Kyu{k,r}'*model.Dinv{k};
                %                 end
                %                 logDetD = logDetD + model.logDetD{k};
                %                 traceDinvyy = traceDinvyy + sum((model.Dinv{k}*model.m{k}).*model.m{k});
                %%%%%%%%%%%%
                model.Dinv{k} = 1./model.D{k};
                model.logDetD{k} = sum(log(model.D{k}));
                for r=1:model.nlf
                    tempo = model.Dinv{k}(:, ones(1, model.k(r)));
                    model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';
                    %                     model.KuyDinv2{r,k} = model.Kyu{k,r}'*sparseDiag(model.Dinv{k});
                end
                logDetD = logDetD + model.logDetD{k};
                traceDinvyy = traceDinvyy + sum((model.Dinv{k}.*model.m{k}).*model.m{k});
            end
        case 'pitc'
            for k =1:model.nout
                KyuKuuinvKuy = zeros(size(model.Kyy{k},1),size(model.Kyy{k},1));
                for r =1: model.nlf,
                    KyuKuuinvKuy = KyuKuuinvKuy + model.Kyu{k,r}*model.KuuinvKuy{r,k};
                end
                model.D{k} = model.Kyy{k} - KyuKuuinvKuy; % In pitc D is a block-diagonal matrix
                %model.D{k} = model.Kyy{k};
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
                        case 1
                            if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                                model.D{k} =(model.D{k} + diag((1/model.beta(k))*(1./model.nRepeats{k})));
                            else
                                error('Model does not contain nRepeats')
                            end
                        case 2
                            % Be aware: the interpretation of beta for this option is as
                            % variance, not precision
                            model.D{k} = (model.D{k} + diag(1./model.beta{k}));
                    end
                else
                    model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
                end
                model.D{k} = checkKernelSymmetry(model.D{k});
                [model.Dinv{k}, model.sqrtD{k}] = pdinv(model.D{k});
                model.logDetD{k} = logdet(model.D{k}, model.sqrtD{k});
                for r=1:model.nlf
                    model.KuyDinv{r,k} = model.Kyu{k,r}'*model.Dinv{k};
                end
                logDetD = logDetD + model.logDetD{k};
                traceDinvyy = traceDinvyy + sum((model.Dinv{k}*model.m{k}).*model.m{k});
            end
            
        case 'dtcvar'
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0                        
                        for k =1:model.nout
                            model.D{k} = 1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
                            model.Dinv{k} = model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
                            model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
                            for r=1:model.nlf
                                model.KuyDinv{r,k} = model.beta(k)*model.Kyu{k,r}';
                            end
                            KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
                            for r =1: model.nlf,
                                KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2);
                            end
                            model.Ktilde{k} = model.Kyy{k} - KyuKuuinvKuy;
                            temp = model.beta(k)*sum(model.Ktilde{k});
                            model.KtildeT = model.KtildeT + temp;
                            logDetD = logDetD + model.logDetD{k};
                            traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                        end
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            for k =1:model.nout
                                model.D{k} = (1/model.beta(k))*(1./model.nRepeats{k});
                                model.Dinv{k} = model.beta(k)*model.nRepeats{k};
                                model.logDetD{k} = sum(log((1/model.beta(k))*(1./model.nRepeats{k})));
                                for r=1:model.nlf
                                    tempo = model.Dinv{k}(:, ones(1, model.k(r)));
                                    model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';                                   
                                end                                
                                KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
                                for r =1: model.nlf,
                                    KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2);
                                end
                                model.Ktilde{k} = model.Kyy{k} - KyuKuuinvKuy;
                                temp = sum(model.beta(k)*model.nRepeats{k}.*model.Ktilde{k});
                                model.KtildeT = model.KtildeT + temp;
                                
                                logDetD = logDetD + model.logDetD{k};
                                traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                            end
                        else
                            error('Model does not contain nRepeats')
                        end
                    case 2
                        for k =1:model.nout
                            % Be aware: the interpretation of beta for this option is as
                            % variance, not precision
                            model.D{k} = 1./model.beta{k};
                            model.Dinv{k} = model.beta{k};
                            model.logDetD{k} = log(1./model.beta{k});
                            for r=1:model.nlf
                                tempo = model.beta{k}(:, ones(1, model.k(r)));
                                model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';                   
                            end                            
                            KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
                            for r =1: model.nlf,
                                KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2);
                            end
                            model.Ktilde{k} = model.Kyy{k} - KyuKuuinvKuy;
                            temp = sum((model.beta{k}).*model.Ktilde{k});
                            model.KtildeT = model.KtildeT + temp;                            
                            logDetD = logDetD + model.logDetD{k};
                            traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
                        end
                end
            else
                for k =1:model.nout
                    model.D{k} = 1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
                    model.Dinv{k} = model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
                    model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
                    for r=1:model.nlf
                        model.KuyDinv{r,k} = model.beta(k)*model.Kyu{k,r}';
                    end
                    KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
                    for r =1: model.nlf,
                        KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2);
                    end
                    model.Ktilde{k} = model.Kyy{k} - KyuKuuinvKuy;
                    temp = model.beta(k)*sum(model.Ktilde{k});
                    model.KtildeT = model.KtildeT + temp;
                    logDetD = logDetD + model.logDetD{k};
                    traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});                   
                end
            end                                   
%             for k =1:model.nout
%                 if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                     switch model.noiseOpt
%                         case 0
%                             model.D{k} = 1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
%                             model.Dinv{k} = model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
%                             model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
%                             for r=1:model.nlf
%                                 model.KuyDinv{r,k} = model.beta(k)*model.Kyu{k,r}';
%                             end
%                         case 1
%                             if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                                 model.D{k} = (1/model.beta(k))*(1./model.nRepeats{k});
%                                 model.Dinv{k} = model.beta(k)*model.nRepeats{k};
%                                 model.logDetD{k} = sum(log((1/model.beta(k))*(1./model.nRepeats{k})));
%                                 for r=1:model.nlf
%                                     tempo = model.Dinv{k}(:, ones(1, model.k(r)));
%                                     model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';
%                                     %model.KuyDinv{r,k} = (diag(model.Dinv{k})*model.Kyu{k,r})';
%                                 end
%                             else
%                                 error('Model does not contain nRepeats')
%                             end
%                         case 2
%                             % Be aware: the interpretation of beta for this option is as
%                             % variance, not precision
%                             model.D{k} = 1./model.beta{k};
%                             model.Dinv{k} = model.beta{k};
%                             model.logDetD{k} = log(1./model.beta{k});
%                             for r=1:model.nlf
%                                 tempo = model.beta{k}(:, ones(1, model.k(r)));
%                                 model.KuyDinv{r,k} = (model.Kyu{k,r}.*tempo)';
%                                 %model.KuyDinv{r,k} = (diag(model.beta{k})*model.Kyu{k,r})';
%                             end
%                     end
%                 else
%                     model.D{k} = 1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
%                     model.Dinv{k} = model.beta(k)*ones(size(model.X{k+model.nlf},1),1);
%                     model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
%                     for r=1:model.nlf
%                         model.KuyDinv{r,k} = model.beta(k)*model.Kyu{k,r}';
%                     end
%                 end
%                 if strcmp(model.approx, 'dtcvar')
%                     KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
%                     for r =1: model.nlf,
%                         KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}.*model.KuuinvKuy{r,k}',2);
%                     end
%                     model.Ktilde{k} = model.Kyy{k} - KyuKuuinvKuy;
%                     if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                         switch model.noiseOpt
%                             case 0
%                                 temp = model.beta(k)*sum(model.Ktilde{k});
%                             case 1
%                                 temp = sum(model.beta(k)*model.nRepeats{k}.*model.Ktilde{k});
%                             case 2
%                                 temp = sum((model.beta{k}).*model.Ktilde{k});
%                         end
%                     else
%                         temp = model.beta(k)*sum(model.Ktilde{k});
%                     end
%                     model.KtildeT = model.KtildeT + temp;
%                 end
%                 logDetD = logDetD + model.logDetD{k};
%                 traceDinvyy = traceDinvyy + sum(model.Dinv{k}.*model.m{k}.*model.m{k});
%             end                        
    end
    
    model.logDetDT = logDetD;
    model.traceDinvyy = traceDinvyy;
    
    for r =1:model.nlf,
        model.KuyDinvy{r,1} = zeros(model.k(r),1);
        for q =1:model.nout,
            model.KuyDinvy{r} = model.KuyDinvy{r} + model.KuyDinv{r,q}*model.m{q};
        end
    end
    
    for k =1:model.nlf,
        KuyDinvKyu = zeros(model.k(k));
        for q =1:model.nout,
            KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}*model.Kyu{q,k};
        end
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            model.A{k,k} = model.KuuGamma{k} + KuyDinvKyu;
        else
            model.A{k,k} = model.Kuu{k} + KuyDinvKyu;
        end
        for r =1:k-1,
            KuyDinvKyu = zeros(model.k(k), model.k(r));
            for q =1:model.nout,
                KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}*model.Kyu{q,r};
            end
            model.A{k,r} = KuyDinvKyu;
            model.A{r,k} = KuyDinvKyu';
        end
    end
    
    A = cell2mat(model.A);
    KuyDinvy = cell2mat(model.KuyDinvy);
    [Ainv, sqrtA]  = pdinv(A);
    sqrtAinv = jitChol(Ainv);
    model.logDetA = logdet(A, sqrtA);
    
    model.Ainv  = mat2cell(Ainv, model.k, model.k);
    model.sqrtA = mat2cell(sqrtA,model.k, model.k);
    model.sqrtAinvKuyDinvy = mat2cell(sqrtAinv*KuyDinvy, model.k, 1);
    
    for r = 1:model.nlf,
        model.AinvKuyDinvy{r,1} = zeros(model.k(r),1);
        for k = 1:model.nlf,
            model.AinvKuyDinvy{r} = model.AinvKuyDinvy{r} + model.Ainv{r,k}*model.KuyDinvy{k};
        end
    end
end


%/~ Old code
% for k =1: model.nout,
%     switch model.approx
%         case {'dtc', 'dtcvar'}
%             if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                 switch model.noiseOpt
%                     case 0
%                         model.D{k} = sparseDiag(1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1));
%                         model.Dinv{k} = sparseDiag(model.beta(k)*ones(size(model.X{k+model.nlf},1),1));
%                         model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
%                     case 1
%                         if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                             model.D{k} = sparseDiag( (1/model.beta(k))*(1./model.nRepeats{k}));
%                             model.Dinv{k} = sparseDiag(model.beta(k)*model.nRepeats{k});
%                             model.logDetD{k} = sum(log((1/model.beta(k))*(1./model.nRepeats{k})));
%                         else
%                             error('Model does not contain nRepeats')
%                         end
%                     case 2
%                         error('Not implemented yet')
%                     case 3
%                         % Be aware: the interpretation of beta for this option is as
%                         % variance, not precision
%                         model.D{k} = sparseDiag(1./model.beta{k});
%                         model.Dinv{k} = sparseDiag(model.beta{k});
%                         model.logDetD{k} = sum(log(1./model.beta{k}));
%                 end
%             else
%                 model.D{k} = sparseDiag(1/model.beta(k)*ones(size(model.X{k+model.nlf},1),1));
%                 model.Dinv{k} = sparseDiag(model.beta(k)*ones(size(model.X{k+model.nlf},1),1));
%                 model.logDetD{k} = -size(model.X{k+model.nlf},1)*log(model.beta(k));
%             end
%             if strcmp(model.approx, 'dtcvar')
%                 KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
%                 for r =1: model.nlf,
%                     KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}*model.KuuinvKuy{r,k}',2);
%                 end
%                 model.Ktilde{k} = model.Kyy{k} - diag(KyuKuuinvKuy);
%                 if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                     switch model.noiseOpt
%                         case 0
%                             temp = model.beta(k)*sum(model.Ktilde{k});
%                         case 1
%                             if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                                 temp = sum(model.beta(k)*model.nRepeats{k}.*model.Ktilde{k});
%                             end
%                         case 2
%
%                         case 3
%                             temp = sum((model.beta{k}).*model.Ktilde{k});
%                     end
%                 else
%                     temp = model.beta(k)*sum(model.Ktilde{k});
%                 end
%                 model.KtildeT = model.KtildeT + temp;
%             end
%         case 'fitc'
%             KyuKuuinvKuy = zeros(size(model.Kyy{k},1),1);
%             for r =1: model.nlf,
%                 KyuKuuinvKuy = KyuKuuinvKuy + sum(model.Kyu{k,r}*model.KuuinvKuy{r,k}',2);
%             end
%             model.D{k} = model.Kyy{k} - diag(KyuKuuinvKuy); % In fitc D is a diagonal matrix.
%             if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                 switch model.noiseOpt
%                     case 0
%                         model.D{k} = (model.D{k} + 1/model.beta(k));
%                     case 1
%                         if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                             model.D{k} = (model.D{k} + (1/model.beta(k))*(1./model.nRepeats{k}));
%                         else
%                             error('Model does not contain nRepeats')
%                         end
%                     case 2
%                         error('Not implemented yet')
%                     case 3
%                         % Be aware: the interpretation of beta for this option is as
%                         % variance, not precision
%                         model.D{k} = (model.D{k} + 1./model.beta{k});
%                 end
%             else
%                 model.D{k} = (model.D{k} + 1/model.beta(k));
%             end
%             model.Dinv{k} = sparseDiag(1./model.D{k});
%             model.logDetD{k} = sum(log(model.D{k}));
%             model.D{k} = sparseDiag(model.D{k});% This is to keep the flow of the gradients
%         case 'pitc'
%             KyuKuuinvKuy = zeros(size(model.Kyy{k},1),size(model.Kyy{k},1));
%             for r =1: model.nlf,
%                 KyuKuuinvKuy = KyuKuuinvKuy + model.Kyu{k,r}*model.KuuinvKuy{r,k};
%             end
%             model.D{k} = model.Kyy{k} - KyuKuuinvKuy; % In pitc D is a block-diagonal matrix
%             %model.D{k} = model.Kyy{k};
%             if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                 switch model.noiseOpt
%                     case 0
%                         model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
%                     case 1
%                         if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                             model.D{k} =(model.D{k} + diag((1/model.beta(k))*(1./model.nRepeats{k})));
%                         else
%                             error('Model does not contain nRepeats')
%                         end
%                     case 2
%                         error('Not implemented yet')
%                     case 3
%                         % Be aware: the interpretation of beta for this option is as
%                         % variance, not precision
%                         model.D{k} = (model.D{k} + diag(1./model.beta{k}));
%                 end
%             else
%                 model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
%             end
%             model.D{k} = checkKernelSymmetry(model.D{k});
%             [model.Dinv{k}, model.sqrtD{k}] = pdinv(model.D{k});
%             model.logDetD{k} = logdet(model.D{k}, model.sqrtD{k});
%         otherwise
%             error('Unknown approximation type')
%     end
%     logDetD = logDetD + model.logDetD{k};
%     traceDinvyy = traceDinvyy + sum((model.Dinv{k}*model.m{k}).*model.m{k});
% end
% model.logDetDT = logDetD;
% model.traceDinvyy = traceDinvyy;
%
% for k =1:model.nlf,
%     for q =1:model.nout,
%         switch model.approx
%             case {'dtc', 'dtcvar'}
%                 if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                     switch model.noiseOpt
%                         case 0
%                             model.KuyDinv{k,q} = model.beta(q)*model.Kyu{q,k}';
%                         case 1
%                             if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                                 model.KuyDinv{k,q} = (diag(model.beta(q)*model.nRepeats{q})*model.Kyu{q,k})';
%                             else
%                                 error('Model does not contain nRepeats')
%                             end
%                         case 2
%                             error('Not implemented yet')
%                         case 3
%                             model.KuyDinv{k,q} = (diag(model.beta{q})*model.Kyu{q,k})';
%                     end
%
%                 else
%                     model.KuyDinv{k,q} = model.beta(q)*model.Kyu{q,k}';
%                 end
%             case {'fitc','pitc'}
%                 model.KuyDinv{k,q} = model.Kyu{q,k}'*model.Dinv{q};
%         end
%     end
% end
%
% for r =1:model.nlf,
%     model.KuyDinvy{r,1} = zeros(model.k(r),1);
%     for q =1:model.nout,
%         model.KuyDinvy{r} = model.KuyDinvy{r} + model.KuyDinv{r,q}*model.m{q};
%     end
% end
%
% for k =1:model.nlf,
%     KuyDinvKyu = zeros(model.k(k));
%     for q =1:model.nout,
%         KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}*model.Kyu{q,k};
%     end
%     if isfield(model, 'gamma') && ~isempty(model.gamma)
%         model.A{k,k} = model.KuuGamma{k} + KuyDinvKyu;
%     else
%         model.A{k,k} = model.Kuu{k} + KuyDinvKyu;
%     end
%     for r =1:k-1,
%         KuyDinvKyu = zeros(model.k(k), model.k(r));
%         for q =1:model.nout,
%             KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}*model.Kyu{q,r};
%         end
%         model.A{k,r} = KuyDinvKyu;
%         model.A{r,k} = KuyDinvKyu';
%     end
% end
%
% A = cell2mat(model.A);
% KuyDinvy = cell2mat(model.KuyDinvy);
% [Ainv, sqrtA]  = pdinv(A);
% sqrtAinv = jitChol(Ainv);
% model.logDetA = logdet(A, sqrtA);
%
% model.Ainv  = mat2cell(Ainv, model.k, model.k);
% model.sqrtA = mat2cell(sqrtA,model.k, model.k);
% model.sqrtAinvKuyDinvy = mat2cell(sqrtAinv*KuyDinvy, model.k, 1);
%
% for r = 1:model.nlf,
%     model.AinvKuyDinvy{r,1} = zeros(model.k(r),1);
%     for k = 1:model.nlf,
%         model.AinvKuyDinvy{r} = model.AinvKuyDinvy{r} + model.Ainv{r,k}*model.KuyDinvy{k};
%     end
% end
% ~/




