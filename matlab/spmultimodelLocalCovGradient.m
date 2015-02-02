function [dLdKyy, dLdKuy, dLdKuu, dLdmu, dLdbeta] = spmultimodelLocalCovGradient(model)

% SPMULTIMODELLOCALCOVGRADIENT Computes the derivatives of the likelihood
% objective function corresponding to the sparse approximation, with
% respect to the kernels of the sparse multi model.
%
% FORMAT
% DESC Computes the derivatives of the likelihood objective function corresponding
% to the sparse approximation
% ARG model : the model structure containing the information about
% the model.
% RETURN dLdKyy: the derivatives of the likelihood with respect to the
% kernels Kyy.
% RETURN dLdKuy: the derivatives of the likelihood with respect to the
% kernels Kuy.
% RETURN dLdKuu: the derivatives of the likelihood with respect to the
% kernels Kuu.
% RETURN dLddmu: the derivatives of the likelihood with respect to the
% mean function.
% RETURN dLdbeta: the derivatives of the likelihood with respect to the
% noise beta parameters.

% COPYRIGHT : Mauricio A Alvarez, 2009

% MULTIGP

Ainv = cell2mat(model.Ainv);

% /~MAURICIO : Last code
% C = mat2cell(Ainv + (AinvKuyDinvy*AinvKuyDinvy'), model.k*ones(1,model.nlf), model.k*ones(1,model.nlf));
% CKuy = mat2cell(cell2mat(C)*Kyu', model.k*ones(1,model.nlf), cellfun('length',model.m));
% ~/

if model.varS
    AinvKuyHatDinvy = cell2mat(model.AinvKuyHatDinvy);
    C = mat2cell(Ainv + (AinvKuyHatDinvy*AinvKuyHatDinvy'), model.k*ones(1,model.nlf),...
         model.k*ones(1,model.nlf));
    dLdKuy = cell(model.nlf,1);
    for r =1:model.nlf,
        tC = cell2mat(C(r,:));
        for k= 1:model.nout,
            exSqSq = model.qs.Sigma(:,r,k) + (model.qs.mean(k,:)')*model.qs.mean(k,r);
            exSqSqBeta = exSqSq*model.beta(k);
            exSqSqBeta2 = kron(exSqSqBeta, ones(model.k, model.sizeX(k)));
            Kuy = cell2mat(model.Kyu(k,:))';
            dLdKuy{r,k} = model.AinvKuyHatDinvy{r}*model.m{k}'*model.beta(k)*model.qs.mean(k,r) + ...
                    - tC*(exSqSqBeta2.*Kuy) + ...
                    model.KuuinvKuy{r, k}*model.beta(k)*model.exS2(k,r);    
        end
    end
    
    dLdKuu = cell(model.nlf,1);
    for r =1:model.nlf,
        KuuinvKuyQKyuKuuinv = zeros(model.k);
        dLdKuu{r} = zeros(model.k);
        for k=1: model.nout,
            KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + ...
                model.KuuinvKuy{r,k}*model.KuuinvKuy{r,k}'*model.beta(k)*model.exS2(k,r);
        end
        dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
    end
    dLdbeta = zeros(1,model.nout);
    dLdmu = cell(model.nout,1);
    dLdKyy = cell(model.nout,model.nlf);
    for k =1:model.nout,        
        temp = kron(-0.5*model.beta(k)*model.exS2(k,:), ...
            sparseDiag(ones(model.sizeX(k),1)));
        dLdKyy(k,:) = mat2cell(temp, model.sizeX(k), ...
            model.sizeX(k)*ones(1, model.nlf));
        Kyu = cell2mat(model.Kyu(k,:));
        exSq = model.qs.mean(k,:)';
        exSq2 = kron(exSq, ones(model.k,1));
        exSqSq = model.qs.Sigma(:,:,k) + (model.qs.mean(k,:)')*(model.qs.mean(k,:));
        exSqSq2 = kron(exSqSq, ones(model.k));
        temp = Kyu*(AinvKuyHatDinvy.*exSq2);
        betaCA = 2*model.m{k}'*temp - sum(sum(exSqSq2.*(Kyu'*Kyu).*cell2mat(C)));         
        dLdmu{k} = model.beta(k)*(model.m{k} - temp);
        dLdbeta(k) = -0.5*(model.Ktilde{k} - (model.sizeX(k)/model.beta(k) ...
            - model.m{k}'*model.m{k} + betaCA));
    end
else
    AinvKuyDinvy = cell2mat(model.AinvKuyDinvy);
    C = mat2cell(Ainv + (AinvKuyDinvy*AinvKuyDinvy'), model.k*ones(1,model.nlf),...
       model.k*ones(1,model.nlf));
    CKuy = cell(model.nlf, model.nout);
    for i =1:model.nlf
        for j =1:model.nout
            whichLatent = model.indLatGivenOut{j};
            CKuy{i,j} = zeros(model.k, model.sizeX(j));
            for k =1:length(whichLatent)
                CKuy{i,j} = CKuy{i,j} +  C{i, whichLatent(k)}*(model.Kyu{j}{k}');
            end
        end
    end
    switch model.approx
        case 'dtcvar'
            dLdKyy = cell(model.nout,1);
            for k=1:model.nout,
                dLdKyy{k} =  -0.5*model.Dinv{k};
            end
            dLdKuy = cell(model.nlf,1);
            for r =1:model.nlf,
                whichOutput = model.indOutGivenLat{r};
                for k= 1:model.numOutGivenLat(r),
                    dLdKuy{r}{k} = model.KuuinvKuy{r}{k}*model.beta(whichOutput(k))...
                        - CKuy{r,whichOutput(k)}*model.beta(whichOutput(k)) + ...
                        model.AinvKuyDinvy{r}*model.m{whichOutput(k)}'*model.beta(whichOutput(k));
                end
            end
            dLdKuu = cell(model.nlf,1);
            for r =1:model.nlf,
                KuuinvKuyQKyuKuuinv = zeros(model.k);
                dLdKuu{r} = zeros(model.k);
                whichOutput = model.indOutGivenLat{r};
                for k=1: model.numOutGivenLat(r),
                    KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + ...
                        model.KuuinvKuy{r}{k}*model.Dinv{whichOutput(k)}*model.KuuinvKuy{r}{k}';
                end
                dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
            end
            dLdbeta = zeros(1,model.nout);
            dLdmu = cell(model.nout,1);
            for k =1:model.nout,
                DinvKyuAinvKuyDinvy = zeros(model.sizeX(k),1);
                H = zeros(model.sizeX(k));
                whichLatent = model.indLatGivenOut{k};
                for r =1:model.numLatGivenOut(k),
                    H = H + model.Kyu{k}{r}*model.AinvKuyDinvy{whichLatent(r)}*model.m{k}'+ ...
                        (model.Kyu{k}{r}*model.AinvKuyDinvy{whichLatent(r)}*model.m{k}')'...
                        - model.Kyu{k}{r}*CKuy{whichLatent(r),k};
                    DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                        model.Dinv{k}*model.Kyu{k}{r}*model.AinvKuyDinvy{whichLatent(r)};
                end
                if ~strcmp(model.kernType,'gg')
                    dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                end
                dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (model.sizeX(k)/model.beta(k) ...
                    - model.m{k}'*model.m{k} + trace(H)));
            end
        case {'fitc','pitc'}
            Q = cell(model.nout,1);
            dLdKyy = cell(model.nout,1);
            dLdmu = cell(model.nout,1);
            for k =1:model.nout,
                DinvKyuAinvKuyDinvy = zeros(model.sizeX(k),1);
                Q{k} = zeros(model.sizeX(k));
                whichLatent = model.indLatGivenOut{k};
                for r =1:model.numLatGivenOut(k),
                    Q{k} = Q{k} + model.Kyu{k}{r}*model.AinvKuyDinvy{whichLatent(r)}*model.m{k}'...
                        + (model.Kyu{k}{r}*model.AinvKuyDinvy{whichLatent(r)}*model.m{k}')' - model.Kyu{k}{r}*CKuy{whichLatent(r),k};
                    DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                        model.Dinv{k}*model.Kyu{k}{r}*model.AinvKuyDinvy{whichLatent(r)};
                end
                Q{k} = model.Dinv{k}*(model.D{k} - model.m{k}*model.m{k}' + Q{k})*model.Dinv{k};
                if strcmp(model.approx, 'fitc'),
                    Q{k} = sparseDiag(diag(Q{k}));
                end
                dLdKyy{k} = -0.5*Q{k};
                if ~strcmp(model.kernType,'gg')
                    dLdmu{k} = model.Dinv{k}*model.m{k} - DinvKyuAinvKuyDinvy;
                end
            end
            dLdKuy = cell(model.nlf,1);
            for r =1:model.nlf,
                whichOutput = model.indOutGivenLat{r};
                for k= 1:model.numOutGivenLat(r),
                    dLdKuy{r}{k} = model.KuuinvKuy{r}{k}*Q{whichOutput(k)}...
                        - CKuy{r,whichOutput(k)}*model.Dinv{whichOutput(k)} + ...
                        model.AinvKuyDinvy{r}*model.m{whichOutput(k)}'*model.Dinv{whichOutput(k)};
                end
            end
            dLdKuu = cell(model.nlf,1);
            for r =1:model.nlf,
                KuuinvKuyQKyuKuuinv = zeros(model.k);
                dLdKuu{r} = zeros(model.k);
                whichOutput = model.indOutGivenLat{r};
                for k=1: model.numOutGivenLat(r),
                    KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + model.KuuinvKuy{r}{k}*Q{whichOutput(k)}*model.KuuinvKuy{r}{k}';
                end
                dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
            end
            if nargout>4
                dLdbeta = zeros(1,model.nout);
                for k =1:model.nout,
                    dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(Q{k});
                end
            end
    end
end
