function model = spmultimodelUpdateAD(model)

% SPMULTIMODELUPDATEAD Update the representations of A and D associated with 
% the sparse multi model.
% FORMAT
% DESC updates the representations of A and D in the model when
% called by spmultimodelExpandParam.
% ARG model : the model for which the representations are being
% updated.
% RETURN model : the model with the A and D representations
% updated.
%
% SEEALSO : spmultimodelExpandParam
%
% COPYRIGHT : Mauricio Alvarez 2009
%
% MODIFICATIONS : Mauricio Alvarez 2010.

% MULTIGP

if model.varS
    for r=1:model.nlf,
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            [model.Kuuinv{r}, model.sqrtKuu{r}, jitter] = pdinv(model.KuuGamma{r});
            model.logDetKuu{r} = logdet(model.KuuGamma{r}, model.sqrtKuu{r});
        else
            [model.Kuuinv{r}, model.sqrtKuu{r}, jitter] = pdinv(model.Kuu{r});
            model.logDetKuu{r} = logdet(model.Kuu{r}, model.sqrtKuu{r});
        end
    end
    for r =1: model.nlf,
        for k =1: model.nout,
            model.KuuinvKuy{r,k} = model.Kuuinv{r}*model.Kyu{k,r}';            
        end
    end
    model.sumCdq = 0;
    model.sumGdqExS2 = 0;
    model.entropyS = 0;
    for k =1: model.nout,
        KyGivenu = 0;
        for r =1: model.nlf,
            tr1 = sum(model.Kyy{k,r});
            tr2 = sum(sum((model.Kyu{k,r}.*(model.KuuinvKuy{r,k})')));
            model.Ed(k,r) = tr1- tr2;
            exS2 = model.qs.Sigma(r,r,k) + (model.qs.mean(k,r))^2;
            KyGivenu = KyGivenu + model.Ed(k,r)*exS2;
            model.sumGdqExS2 = model.sumGdqExS2 + ...
                model.gammas(k+(r-1)*model.nout)*exS2;
            model.exS2(k,r) = exS2;
            
        end
        model.entropyS = model.entropyS + ...
                logdet(model.qs.Sigma(:,:,k)) + model.nlf*(1+log(2*pi)); 
        model.D{k} = sparseDiag(1/model.beta(k)*ones(model.sizeX(k),1));
        model.Dinv{k} = sparseDiag(model.beta(k)*ones(model.sizeX(k),1));
        model.logDetD{k} = - model.sizeX(k)*log(model.beta(k));
        model.Ktilde{k} = KyGivenu;
        model.sumCdq = model.sumCdq + model.beta(k)*KyGivenu;
    end

    for r =1:model.nlf,
        model.KuyHatDinvy{r,1} = zeros(model.k,1);        
        for q =1:model.nout,
            model.KuyHatDinvy{r} = model.KuyHatDinvy{r} + ...
                model.beta(q)*model.qs.mean(q,r)*((model.Kyu{q,r})'*model.m{q});
        end
    end

    % Compute A
    for k =1:model.nlf,
        Kuy1 = (cell2mat(model.Kyu(:,k)))';
        exS = squeeze(model.qs.Sigma(k,k,:)) + model.qs.mean(:,k).^2;
        exS = exS.*model.beta'; % Times beta
        exS2 = kron(exS, ones(model.sizeX(1),1));
        exSM = repmat(exS2', model.k, 1);
        Kuy1exSM = Kuy1.*exSM;
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            model.A{k,k} = model.KuuGamma{k} + Kuy1exSM*(Kuy1');
        else
            model.A{k,k} = model.Kuu{k} + Kuy1exSM*(Kuy1');
        end       
        for r =1:k-1,
            Kyu2 = cell2mat(model.Kyu(:,r));
            exS = squeeze(model.qs.Sigma(k,r,:)) + ...
                model.qs.mean(:,k).*model.qs.mean(:,r);
            exS = exS.*model.beta'; % Times beta
            exS2 = kron(exS, ones(model.sizeX(1),1));
            exSM = repmat(exS2', model.k, 1);
            Kuy1exSM = Kuy1.*exSM;
            model.A{k,r} = Kuy1exSM*Kyu2;
            model.A{r,k} = (model.A{k,r})';
        end
    end

    A = sparse(cell2mat(model.A));
    %Ainv = real(inv(A));
    %model.logDetA = real(log(det(A)));
    Ainv = pdinv(A); 
    model.logDetA = logdet(A);
    model.Ainv  = mat2cell(Ainv, model.k*ones(1, model.nlf), model.k*ones(1, model.nlf));
    
    for r = 1:model.nlf,
        model.AinvKuyHatDinvy{r,1} = zeros(model.k,1);
        for k = 1:model.nlf,
            model.AinvKuyHatDinvy{r} = model.AinvKuyHatDinvy{r} + model.Ainv{r,k}*model.KuyHatDinvy{k};
        end
    end
    
else

    for r=1:model.nlf,
        % Kuuinv really is the inverse of Kuu + regularization term
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            [model.Kuuinv{r}, model.sqrtKuu{r}, jitter] = pdinv(model.KuuGamma{r});
            model.logDetKuu{r} = logdet(model.KuuGamma{r}, model.sqrtKuu{r});
        else
            [model.Kuuinv{r}, model.sqrtKuu{r}, jitter] = pdinv(model.Kuu{r});
            model.logDetKuu{r} = logdet(model.Kuu{r}, model.sqrtKuu{r});
        end
    end

    for r =1: model.nlf,
        for k =1: model.numOutGivenLat(r),
            %model.KuuinvKuy{r,k} = model.Kuuinv{r}*model.Kyu{k,r}';
            model.KuuinvKuy{r,1}{k,1} = model.Kuuinv{r}*model.Kuy{r}{k};
        end
    end

    for k =1: model.nout,
        KyuKuuinvKuy = zeros(model.sizeX(k));
        %KyuKuuinvKuy2 = zeros(model.sizeX(k));
        whichLatent = model.indLatGivenOut{k};
        for r =1: model.numLatGivenOut(k),
            whichOutputs = model.indOutGivenLat{whichLatent(r)};
            subsetLatent = find(whichOutputs == k);
            KyuKuuinvKuy = KyuKuuinvKuy + model.Kyu{k}{r}*model.KuuinvKuy{whichLatent(r)}{subsetLatent};
            %KyuKuuinvKuy2 = KyuKuuinvKuy2 + model.Kyu{k}{r}*model.Kuuinv{whichLatent(r)}*model.Kyu{k}{r}';
            %KyuKuuinvKuy = KyuKuuinvKuy + model.Kyu{k,r}*model.KuuinvKuy{r,k};
        end
        switch model.approx
            case 'dtcvar'
                model.D{k} = sparseDiag(1/model.beta(k)*ones(model.sizeX(k),1));
                model.Dinv{k} = sparseDiag(model.beta(k)*ones(model.sizeX(k),1));
                model.logDetD{k} = - model.sizeX(k)*log(model.beta(k));
                model.Ktilde{k} = model.Kyy{k} - diag(KyuKuuinvKuy);
            case 'fitc'
                model.D{k} = model.Kyy{k} - diag(KyuKuuinvKuy); % In fitc D is a diagonal matrix.
                model.D{k} = (model.D{k} + 1/model.beta(k));
                model.Dinv{k} = sparseDiag(1./model.D{k});
                model.logDetD{k} = sum(log(model.D{k}));
                model.D{k} = sparseDiag(model.D{k});% This is to keep the flow of the gradients
            case 'pitc'
                model.D{k} = model.Kyy{k} - KyuKuuinvKuy; % In pitc D is a block-diagonal matrix
                %   model.D2{k} = model.Kyy{k} - KyuKuuinvKuy2; % In pitc D is a block-diagonal matrix
                model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
                model.D{k} = checkKernelSymmetry(model.D{k});
                [model.Dinv{k}, model.sqrtD{k}, jitter] = pdinv(model.D{k});
                model.logDetD{k} = logdet(model.D{k}, model.sqrtD{k});
            otherwise
                error('Unknown approximation type')
        end

    end

    for r =1:model.nlf,
        whichOutputs = model.indOutGivenLat{r};
        for q =1:model.numOutGivenLat(r),
            switch model.approx
                case 'dtcvar'
                    model.KuyDinv{r,1}{q,1} = model.Kuy{r}{q}*model.beta(whichOutputs(q));
                case {'fitc','pitc'}
                    model.KuyDinv{r,1}{q,1} = model.Kuy{r}{q}*model.Dinv{whichOutputs(q)};
            end
        end
    end

    for r =1:model.nlf,
        model.KuyDinvy{r,1} = zeros(model.k,1);
        whichOutputs = model.indOutGivenLat{r};
        for q =1:model.numOutGivenLat(r),
            model.KuyDinvy{r} = model.KuyDinvy{r} + model.KuyDinv{r}{q}*model.m{whichOutputs(q)};
        end
    end

    % Compute A
    % First the terms in the main diagonal
    for k =1:model.nlf,
        KuyDinvKyu = zeros(model.k);
        for q =1:model.numOutGivenLat(k),
            KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k}{q}*(model.Kuy{k}{q})';
        end
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            model.A{k,k} = model.KuuGamma{k} + KuyDinvKyu;
        else
            model.A{k,k} = model.Kuu{k} + KuyDinvKyu;
        end
    end

    % Second the terms off the diagonal

    [rowsNonzero, colsNonzero] = find(model.latConnect  ~= 0);

    for i = 1:length(rowsNonzero)
        KuyDinvKyu = zeros(model.k);
        if rowsNonzero(i)~= colsNonzero(i)
            whichLatent1 = model.indOutGivenLat{rowsNonzero(i)};
            whichLatent2 = model.indOutGivenLat{colsNonzero(i)};
            cont = 0;
            subwhichLatent1 = zeros(1,model.latConnect(rowsNonzero(i),colsNonzero(i)));
            subwhichLatent2 = zeros(1,model.latConnect(rowsNonzero(i),colsNonzero(i)));
            if length(whichLatent2)> length(whichLatent1)
                for j=1:length(whichLatent2)
                    index = find(whichLatent2(j) == whichLatent1);
                    if ~isempty(index)
                        cont = cont + 1;
                        subwhichLatent1(cont) =  index;
                        subwhichLatent2(cont) =  j;
                    end
                end
            else
                for j=1:length(whichLatent1)
                    index = find(whichLatent1(j) == whichLatent2);
                    if ~isempty(index)
                        cont = cont + 1;
                        subwhichLatent1(cont) =  j;
                        subwhichLatent2(cont) =  index;
                    end
                end
            end
            if length(subwhichLatent1)~=model.latConnect(rowsNonzero(i),colsNonzero(i))
                error('Fatal error: the number of forces to multiply is different')
            end
            for q =1:length(subwhichLatent1),
                KuyDinvKyu = KuyDinvKyu +  model.KuyDinv{rowsNonzero(i)}{subwhichLatent1(q)}...
                    *(model.Kuy{colsNonzero(i)}{subwhichLatent2(q)})';
            end
            model.A{rowsNonzero(i),colsNonzero(i)} = KuyDinvKyu;
            model.A{colsNonzero(i),rowsNonzero(i)} = KuyDinvKyu';
        end
    end

    % Finally fill the empty ones with zeros

    [rowsZero, colsZero] = find(model.latConnect == 0);

    for i =1:size(rowsZero,1)
        model.A{rowsZero(i),colsZero(i)} = zeros(model.k);
    end

    A = sparse(cell2mat(model.A));
    %Ainv = real(inv(A));
    %model.logDetA = real(log(det(A)));
    Ainv = pdinv(A);
    model.logDetA = logdet(A);
    model.Ainv  = mat2cell(Ainv, model.k*ones(1, model.nlf), model.k*ones(1, model.nlf));

    [rowsNonzero, colsNonzero] = find(model.latConnectInv ~= 0);

    for r = 1:model.nlf,
        model.AinvKuyDinvy{r,1} = zeros(model.k,1);
        % index = find(rowsNonzero == r);
        colsIndex = colsNonzero(rowsNonzero == r);
        for k = 1:length(colsIndex),
            model.AinvKuyDinvy{r} = model.AinvKuyDinvy{r} + model.Ainv{r,colsIndex(k)}*model.KuyDinvy{colsIndex(k)};
        end
    end
end
