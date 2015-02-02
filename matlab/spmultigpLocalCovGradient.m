function [dLdKyy, dLdKuy, dLdKuu, dLdmu, dLdbeta] = spmultigpLocalCovGradient(model)

% SPMULTIGPLOCALCOVGRADIENT Computes the derivatives of the likelihood
% objective function corresponding to the sparse approximation, with
% respect to the kernels of the multi output gp.
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

% COPYRIGHT : Mauricio A Alvarez, 2008, 2009, 2010

% MULTIGP


if strcmp(model.kernType, 'lmc')
    AinvKuyDinvy = cell2mat(model.AinvKuyDinvy2);
    CM = model.Ainv + AinvKuyDinvy*AinvKuyDinvy';
    startOne = 1;
    endOne = 0;

    A = cell(model.nlf, 1);
    C = cell(model.nlf);
    for r=1:model.nlf
        endOne = endOne + model.k(r)*model.rankCorregMatrix;
        subCM = CM(startOne:endOne, startOne:endOne);
        C{r,r} = mat2cell(subCM, model.k(r)*ones(model.rankCorregMatrix,1), ...
            model.k(r)*ones(model.rankCorregMatrix,1));
        startTwo = 1;
        endTwo = 0;
        for k =1:r-1
            endTwo = endTwo + model.k(k)*model.rankCorregMatrix;
            subCM = CM(startOne:endOne, startTwo:endTwo);
            C{r, k} = mat2cell(subCM, model.k(r)*ones(model.rankCorregMatrix,1), ...
                model.k(k)*ones(model.rankCorregMatrix,1));
            subCM = CM(startTwo:endTwo, startOne:endOne);
            C{k, r} = mat2cell(subCM, model.k(k)*ones(model.rankCorregMatrix,1), ...
                model.k(r)*ones(model.rankCorregMatrix,1));
            startTwo = endTwo + 1;
        end
        A{r} = model.kern.comp{r}.comp{model.nlf+1}.A;
        startOne = endOne + 1;
    end
    CKuy = cell(model.nlf*model.rankCorregMatrix, model.nout);
    for r=1:model.nlf
        for i=1:model.rankCorregMatrix
            linI = (r-1)*model.rankCorregMatrix + i;
            for j =1:model.nout
                CKuy{linI, j} = zeros(model.k(r), size(model.X{model.nlf+j},1));
                for n=1:model.nlf
                    for m=1:model.rankCorregMatrix
                        CKuy{linI, j} = CKuy{linI, j} + ...
                            C{r,n}{i,m}*(A{n}(j,m)*model.Kyu{j,n}');
                    end
                end
            end
        end
    end
    switch model.approx
        case 'dtc'
            dLdKyy = cell(model.nout,1);
            dLdKuy = cell(model.nlf*model.rankCorregMatrix,model.nout);
            for r =1:model.nlf,
                for i= 1:model.rankCorregMatrix
                    linI = (r-1)*model.rankCorregMatrix + i;
                    for k= 1:model.nout,
                        if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                                model.noiseOpt == 1
                            if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                                dLdKuy{linI,k} =  - CKuy{linI,k}*diag(model.beta(k)*model.nRepeats{k}) ...
                                    + model.AinvKuyDinvy{r}(:,i)*model.m{k}'*diag(model.beta(k)*model.nRepeats{k});
                            else
                                error('No information about number of repetitions provided')
                            end
                        else
                            dLdKuy{linI,k} = - CKuy{linI,k}*model.beta(k) + ...
                                model.AinvKuyDinvy{r}(:,i)*model.m{k}'*model.beta(k);
                        end
                    end
                end
            end
            dLdKuu = cell(model.nlf*model.rankCorregMatrix,1);
            for r =1:model.nlf,
                for i=1:model.rankCorregMatrix
                    linI = (r-1)*model.rankCorregMatrix + i;
                    dLdKuu{linI} = 0.5*((model.Kuuinv{r} - C{r,r}{i,i}));
                end
            end
            dLdbeta = zeros(1,model.nout);
            dLdmu = cell(model.nout,1);
            for k =1:model.nout,
                DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                H = zeros(size(model.X{k+model.nlf},1));
                for r =1:model.nlf,
                    for i=1:model.rankCorregMatrix
                        linI = (r-1)*model.rankCorregMatrix + i;
                        H = H + A{r}(k,i)*model.Kyu{k,r}*model.AinvKuyDinvy{r}(:,i)*model.m{k}'...
                            + (A{r}(k,i)*model.Kyu{k,r}*model.AinvKuyDinvy{r}(:,i)*model.m{k}')' ...
                            - A{r}(k,i)*model.Kyu{k,r}*CKuy{linI,k};
                        DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                            model.KuyDinv{r,k}(:,:,i)'*model.AinvKuyDinvy{r}(:,i);
                    end
                end
                if ~strcmp(model.kernType,'gg')
                    if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                            model.noiseOpt == 1
                        dLdmu{k} = (model.beta(k)*model.nRepeats{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
                    else
                        dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                    end
                end
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                        model.noiseOpt == 1
                    dLdbeta(k) = -0.5*model.beta(k)^(-2)*trace(sparseDiag((1./model.nRepeats{k}))*...
                        (model.Dinv{k}*( - ...
                        (model.D{k} - model.m{k}*model.m{k}' + H))*model.Dinv{k}));
                else
                    dLdbeta(k) = -0.5*(- (size(model.X{k+model.nlf},1)*1/model.beta(k)...
                        - model.m{k}'*model.m{k} + trace(H)));
                end
            end
        case {'fitc','pitc'}
            Q = cell(model.nout,1);
            dLdKyy = cell(model.nout,1);
            dLdmu = cell(model.nout,1);
            for k =1:model.nout,
                DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                Q{k} = zeros(size(model.X{k+model.nlf},1));
                for r =1:model.nlf,
                    for i=1:model.rankCorregMatrix
                        linI = (r-1)*model.rankCorregMatrix + i;
                        Q{k} = Q{k} + A{r}(k,i)*model.Kyu{k,r}*model.AinvKuyDinvy{r}(:,i)*model.m{k}'...
                            + (A{r}(k,i)*model.Kyu{k,r}*model.AinvKuyDinvy{r}(:,i)*model.m{k}')' ...
                            - A{r}(k,i)*model.Kyu{k,r}*CKuy{linI,k};
                    end
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
            dLdKuy = cell(model.nlf*model.rankCorregMatrix,model.nout);
            for r =1:model.nlf,
                for i= 1:model.rankCorregMatrix
                    linI = (r-1)*model.rankCorregMatrix + i;
                    for k= 1:model.nout,
                        dLdKuy{linI,k} = A{r}(k,i)*model.KuuinvKuy{r,k}*Q{k}...
                            - CKuy{linI,k}*model.Dinv{k} + model.AinvKuyDinvy{r}(:,i)*model.m{k}'*model.Dinv{k};
                    end
                end
            end
            dLdKuu = cell(model.nlf*model.rankCorregMatrix,1);
            for r =1:model.nlf,
                for i=1:model.rankCorregMatrix
                    linI = (r-1)*model.rankCorregMatrix + i;
                    KuuinvKuyQKyuKuuinv = zeros(model.k(r));
                    for k=1: model.nout,
                        KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + ...
                            (A{r}(k,i)^2)*model.KuuinvKuy{r,k}*Q{k}*model.KuuinvKuy{r,k}';
                    end
                    dLdKuu{linI} = 0.5*((model.Kuuinv{r} - C{r,r}{i,i}) - KuuinvKuyQKyuKuuinv);
                end
            end
            if nargout>4
                dLdbeta = zeros(1,model.nout);
                for k =1:model.nout,
                    if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                            (model.noiseOpt == 1)
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(diag((1./model.nRepeats{k}))*Q{k});
                        end
                    else
                        dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(Q{k});
                    end
                end
            end
        case 'dtcvar'
            dLdKyy = cell(model.nout,1);
            for k=1:model.nout,
                dLdKyy{k} =  -0.5*model.Dinv{k};
            end
            dLdKuy = cell(model.nlf*model.rankCorregMatrix,model.nout);
            for r =1:model.nlf,
                for i=1:model.rankCorregMatrix
                    linI = (r-1)*model.rankCorregMatrix + i;
                    for k= 1:model.nout,
                        if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                                model.noiseOpt == 1
                            if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                                dLdKuy{linI,k} = A{r}(k,i)*model.KuuinvKuy{r,k}*diag(model.beta(k)*model.nRepeats{k})...
                                    - CKuy{linI,k}*diag(model.beta(k)*model.nRepeats{k}) ...
                                    + model.AinvKuyDinvy{r}(:,i)*model.m{k}'*diag(model.beta(k)*model.nRepeats{k});
                            else
                                error('No information about number of repetitions provided')
                            end
                        else
                            dLdKuy{linI,k} = A{r}(k,i)*model.KuuinvKuy{r,k}*model.beta(k)...
                                - CKuy{linI,k}*model.beta(k) + model.AinvKuyDinvy{r}(:,i)*model.m{k}'*model.beta(k);
                        end
                    end
                end
            end
            dLdKuu = cell(model.nlf*model.rankCorregMatrix,1);
            for r =1:model.nlf,
                for i=1:model.rankCorregMatrix
                    linI = (r-1)*model.rankCorregMatrix + i;
                    KuuinvKuyQKyuKuuinv = zeros(model.k(r));
                    for k=1: model.nout,
                        KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv +  ...
                            (A{r}(k,i)^2)*model.KuuinvKuy{r,k}*model.Dinv{k}*model.KuuinvKuy{r,k}';
                    end
                    dLdKuu{linI} = 0.5*((model.Kuuinv{r} - C{r,r}{i,i}) - KuuinvKuyQKyuKuuinv);
                end
            end
            dLdbeta = zeros(1,model.nout);
            dLdmu = cell(model.nout,1);
            for k =1:model.nout,
                DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                H = zeros(size(model.X{k+model.nlf},1));
                for r =1:model.nlf,
                    for i=1:model.rankCorregMatrix
                        linI = (r-1)*model.rankCorregMatrix + i;
                        H = H + A{r}(k,i)*model.Kyu{k,r}*model.AinvKuyDinvy{r}(:,i)*model.m{k}' ...
                              + (A{r}(k,i)*model.Kyu{k,r}*model.AinvKuyDinvy{r}(:,i)*model.m{k}')'...
                              - A{r}(k,i)*model.Kyu{k,r}*CKuy{linI,k};
                        DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                            A{r}(k,i)*model.KuyDinv{r,k}(:,:,i)'*model.AinvKuyDinvy{r}(:,i);
                    end
                end
                if ~strcmp(model.kernType,'gg')
                    if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                        model.noiseOpt == 1
                    else
                        dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                    end
                end
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && ...
                        model.noiseOpt == 1
                    dLdbeta(k) = -0.5*model.beta(k)^(-2)*trace(sparseDiag((1./model.nRepeats{k}))*...
                        (model.Dinv{k}*(sparseDiag(model.Ktilde{k}) - ...
                        (model.D{k} - model.m{k}*model.m{k}' + H))*model.Dinv{k}));
                else
                    dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*...
                        1/model.beta(k) - model.m{k}'*model.m{k} + trace(H))  );
                end
            end
    end
else

    Ainv = cell2mat(model.Ainv);
    AinvKuyDinvy = cell2mat(model.AinvKuyDinvy);
    Kyu = cell2mat(model.Kyu);
    C = mat2cell(Ainv + (AinvKuyDinvy*AinvKuyDinvy'), model.k, model.k);
    CKuy = mat2cell(cell2mat(C)*Kyu', model.k, cellfun('length',model.m));

    switch model.approx
        case 'dtc'
            dLdKyy = cell(model.nout,1);
            dLdKuy = cell(model.nlf,model.nout);
            dLdbeta = zeros(1,model.nout);
            dLdmu = cell(model.nout,1);
            dLdKuu = cell(model.nlf,1);
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        for r =1:model.nlf,
                            for k= 1:model.nout,
                                dLdKuy{r,k} = - CKuy{r,k}*model.beta(k) + ...
                                    model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                            end
                        end                        
                        for k =1:model.nout,                            
                            H = zeros(size(model.X{k+model.nlf},1),1);
                            for r =1:model.nlf,
                                temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};
                                H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);                                
                            end
                            dLdbeta(k) = -0.5*( - (size(model.X{k+model.nlf},1)*1/model.beta(k) - ...
                                model.m{k}'*model.m{k} + sum(H)));
                        end
                        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                            for k =1:model.nout,
                                DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                                for r =1:model.nlf,                                
                                    DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                                        model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
                                end
                                dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                            end
                        end                                                    
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            for r =1:model.nlf,
                                for k= 1:model.nout,
                                    tempo1 = (model.beta(k)*model.nRepeats{k})';
                                    tempo2 = tempo1(ones(model.k(r),1), :);
                                    dLdKuy{r,k} =  - CKuy{r,k}.*tempo2 + (model.AinvKuyDinvy{r}*model.m{k}').*tempo2;
                                end
                            end
                            for k =1:model.nout,                                
                                H = zeros(size(model.X{k+model.nlf},1),1);
                                for r =1:model.nlf,
                                    temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};
                                    H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);                                    
                                end
                                dLdbeta(k) = -0.5*model.beta(k)^(-2)*sum((1./model.nRepeats{k}).*...
                                    (model.Dinv{k}.*( - (model.D{k} - model.m{k}.*model.m{k} + H)).*model.Dinv{k}));
                            end                                                                
                        else
                            error('Information about number of repetitions is not provided')
                        end
                    case 2
                        error('Not implemented yet')                        
                end
            else
                for r =1:model.nlf,
                    for k= 1:model.nout,
                         dLdKuy{r,k} = - CKuy{r,k}*model.beta(k) + ...
                            model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                    end
                end
                for k =1:model.nout,
                    H = zeros(size(model.X{k+model.nlf},1),1);
                    for r =1:model.nlf,
                        temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};
                        H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);
                    end
                     dLdbeta(k) = -0.5*(- (size(model.X{k+model.nlf},1)*1/model.beta(k)...
                        - model.m{k}'*model.m{k} + sum(H)));
                end
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                    for k =1:model.nout,
                        DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                        for r =1:model.nlf,
                            DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                                model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
                        end
                        dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                    end
                end
            end           
            for r =1:model.nlf,
                dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}));
            end                        
%             for r =1:model.nlf,
%                 for k= 1:model.nout,
%                     if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                         switch model.noiseOpt
%                             case 0
%                                 dLdKuy{r,k} = - CKuy{r,k}*model.beta(k) + ...
%                                     model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
%                             case 1
%                                 if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                                     tempo1 = (model.beta(k)*model.nRepeats{k})';
%                                     tempo2 = tempo1(ones(model.k(r),1), :);
%                                     dLdKuy{r,k} =  - CKuy{r,k}.*tempo2 + (model.AinvKuyDinvy{r}*model.m{k}').*tempo2;
%                                 else
%                                     error('Information about number of repetitions is not provided')
%                                 end
%                             case 2
%                                 tempo1 = (model.beta{k})';
%                                 tempo2 = tempo1(ones(model.k(r),1), :);
%                                 dLdKuy{r,k} =  - CKuy{r,k}.*tempo2 + (model.AinvKuyDinvy{r}*model.m{k}').*tempo2;
%                         end
%                     else
%                         dLdKuy{r,k} = - CKuy{r,k}*model.beta(k) + ...
%                             model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
%                     end
%                 end
%             end
%             dLdKuu = cell(model.nlf,1);
%             for r =1:model.nlf,
%                 dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}));
%             end
%             dLdbeta = zeros(1,model.nout);
%             dLdmu = cell(model.nout,1);
%             for k =1:model.nout,
%                 DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
%                 H = zeros(size(model.X{k+model.nlf},1),1);
%                 for r =1:model.nlf,
%                     temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};       
%                     H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);                    
%                     DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
%                         model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
%                 end
%                 if ~strcmp(model.kernType,'gg')
%                     if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                         switch model.noiseOpt
%                             case 0
%                                 dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
%                             case 1
%                                 dLdmu{k} = (model.beta(k)*model.nRepeats{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
%                             case 2
%                                 dLdmu{k} = (model.beta{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
%                         end
%                     else
%                         dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
%                     end
%                 end
%                 if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                     switch model.noiseOpt
%                         case 0
%                             dLdbeta(k) = -0.5*( - (size(model.X{k+model.nlf},1)*1/model.beta(k) - ...
%                                 model.m{k}'*model.m{k} + sum(H)));
%                         case 1
%                             dLdbeta(k) = -0.5*model.beta(k)^(-2)*sum((1./model.nRepeats{k}).*...
%                                 (model.Dinv{k}.*( - ...
%                                 (model.D{k} - model.m{k}.*model.m{k} + H)).*model.Dinv{k}));
%                         case 2
%                             error('Not implemented yet')
%                     end
%                 else
%                     dLdbeta(k) = -0.5*(- (size(model.X{k+model.nlf},1)*1/model.beta(k)...
%                         - model.m{k}'*model.m{k} + sum(H)));
%                 end
%             end
        case 'fitc'
            Q = cell(model.nout,1);            
            dLdKyy = cell(model.nout,1);
            dLdmu = cell(model.nout,1);
            Dinvy = cell(model.nout,1);
            for k =1:model.nout,
                Q{k} = zeros(size(model.X{k+model.nlf},1),1);
                for r =1:model.nlf,
                    KyuAinvKuyDinvy = model.Kyu{k,r}*model.AinvKuyDinvy{r};
                    KyuAinvKuyDinvyy = KyuAinvKuyDinvy.*model.m{k};
                    Q{k} = Q{k} + 2*KyuAinvKuyDinvyy - sum(model.Kyu{k,r}.*CKuy{r,k}',2);
                end
                Dinvy{k} = model.Dinv{k}.*model.m{k}; 
                Q{k} = model.Dinv{k}- Dinvy{k}.*Dinvy{k} + model.Dinv{k}.*Q{k}.*model.Dinv{k};
                if isfield(model, 'useKernDiagGradient') && model.useKernDiagGradient
                    dLdKyy{k} = -0.5*Q{k};
                else
                    dLdKyy{k} = -0.5*sparseDiag(Q{k});
                end                
            end
            dLdKuy = cell(model.nlf,model.nout);
            KuuinvKuyQ = cell(model.nlf,model.nout);
            for r =1:model.nlf,
                for k= 1:model.nout,
                    tempoQ = Q{k}(:,ones(1,model.k(r)));
                    tempoD = model.Dinv{k}(:,ones(1,model.k(r)));
                    KuuinvKuyQ{r,k} = model.KuuinvKuy{r,k}.*tempoQ';
                    dLdKuy{r,k} = KuuinvKuyQ{r,k} - CKuy{r,k}.*tempoD'  + model.AinvKuyDinvy{r}*Dinvy{k}';
                    %KuuinvKuyQ{r,k} = model.KuuinvKuy{r,k}*diag(Q{k});
                    %dLdKuy{r,k} = KuuinvKuyQ2{r,k} - CKuy{r,k}*diag(model.Dinv{k}) ...
                    %    + model.AinvKuyDinvy{r}*Dinvy{k}';
                end
            end
            dLdKuu = cell(model.nlf,1);
            for r =1:model.nlf,
                KuuinvKuyQKyuKuuinv = zeros(model.k(r));
                dLdKuu{r} = zeros(model.k(r));
                for k=1: model.nout,
                     KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv ...
                        + KuuinvKuyQ{r,k}*model.KuuinvKuy{r,k}';
                end
                dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
            end
            if nargout>4
                dLdbeta = zeros(1,model.nout);           
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    for k =1:model.nout,
                        switch model.noiseOpt
                            case 0
                                dLdbeta(k) = 0.5*model.beta(k)^(-2)*sum(Q{k});
                            case 1
                                dLdbeta(k) = 0.5*model.beta(k)^(-2)*sum((1./model.nRepeats{k}).*Q{k});
                            case 2
                                error('Not implemented yet')
                        end
                    end
                else
                    for k =1:model.nout,
                        dLdbeta(k) = 0.5*model.beta(k)^(-2)*sum(Q{k});
                    end
                end            
            end        
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                for k =1:model.nout,
                    DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                    for r =1:model.nlf,
                        DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                            model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
                    end
                    dLdmu{k} = Dinvy{k} - DinvKyuAinvKuyDinvy;
                end
            end            
        case 'pitc'
            Q = cell(model.nout,1);
            dLdKyy = cell(model.nout,1);
            dLdmu = cell(model.nout,1);
            Dinvy = cell(model.nout,1);
            for k =1:model.nout,
                Q{k} = zeros(size(model.X{k+model.nlf},1));
                for r =1:model.nlf,
                    KyuAinvKuyDinvy = model.Kyu{k,r}*model.AinvKuyDinvy{r};
                    KyuAinvKuyDinvyy = KyuAinvKuyDinvy*model.m{k}';
                    Q{k} = Q{k} + KyuAinvKuyDinvyy + KyuAinvKuyDinvyy' - model.Kyu{k,r}*CKuy{r,k};
                end
                Dinvy{k} = model.Dinv{k}*model.m{k}; 
                Q{k} = model.Dinv{k}- Dinvy{k}*Dinvy{k}' + model.Dinv{k}*Q{k}*model.Dinv{k};
                dLdKyy{k} = -0.5*Q{k};
            end
            dLdKuy = cell(model.nlf,model.nout);
            KuuinvKuyQ = cell(model.nlf,model.nout);
            for r =1:model.nlf,
                for k= 1:model.nout,
                    KuuinvKuyQ{r,k} = model.KuuinvKuy{r,k}*Q{k};
                    dLdKuy{r,k} = KuuinvKuyQ{r,k} - CKuy{r,k}*model.Dinv{k} ...
                        + model.AinvKuyDinvy{r}*Dinvy{k}';
                end
            end
            dLdKuu = cell(model.nlf,1);
            for r =1:model.nlf,
                KuuinvKuyQKyuKuuinv = zeros(model.k(r));
                dLdKuu{r} = zeros(model.k(r));
                for k=1: model.nout,
                     KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv ...
                        + KuuinvKuyQ{r,k}*model.KuuinvKuy{r,k}';
                end
                dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
            end
            if nargout>4
                dLdbeta = zeros(1,model.nout);
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    for k =1:model.nout,
                        switch model.noiseOpt
                            case 0
                                dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(Q{k});
                            case 1
                                dLdbeta(k) = 0.5*model.beta(k)^(-2)*...
                                    trace(diag((1./model.nRepeats{k}))*Q{k});
                            case 2
                                error('Not implemented yet')
                        end
                    end
                else
                    for k =1:model.nout,
                        dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(Q{k});
                    end
                end 
            end
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                for k =1:model.nout,
                    DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                    for r =1:model.nlf,
                        DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                            model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
                    end
                    dLdmu{k} = Dinvy{k} - DinvKyuAinvKuyDinvy;
                end
            end   
            
        case 'dtcvar'            
            dLdKyy = cell(model.nout,1);
            if isfield(model, 'useKernDiagGradient') && model.useKernDiagGradient
                for k=1:model.nout,
                    dLdKyy{k} =  -0.5*model.Dinv{k};
                end
            else
                for k=1:model.nout,
                    dLdKyy{k} =  -0.5*sparseDiag(model.Dinv{k});
                end
            end
            dLdKuy = cell(model.nlf,model.nout);
            dLdbeta = zeros(1,model.nout);
            dLdmu = cell(model.nout,1);
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        for r =1:model.nlf,
                            for k= 1:model.nout,
                                dLdKuy{r,k} = model.KuuinvKuy{r,k}*model.beta(k)...
                                    - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                                
                            end
                        end
                        for k =1:model.nout,
                            H = zeros(size(model.X{k+model.nlf},1),1);
                            for r =1:model.nlf,
                                temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};
                                H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);                                
                            end                            
                            dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*1/model.beta(k) - ...
                                model.m{k}'*model.m{k} + sum(H)));                            
                        end                        
                        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)                    
                            for k =1:model.nout,
                                DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                                for r =1:model.nlf,
                                    DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                                        model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
                                end
                                dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                            end
                        end      
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            for r =1:model.nlf,
                                for k= 1:model.nout,
                                    tempo1 = (model.beta(k)*model.nRepeats{k})';
                                    tempo2 = tempo1(ones(model.k(r),1), :);
                                    dLdKuy{r,k} = (model.KuuinvKuy{r,k} - CKuy{r,k} + (model.AinvKuyDinvy{r}*model.m{k}')).*tempo2;
                                end
                            end
                            for k =1:model.nout,                                
                                H = zeros(size(model.X{k+model.nlf},1),1);
                                for r =1:model.nlf,
                                    temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};
                                    H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);                                    
                                end
                                dLdbeta(k) = -0.5*model.beta(k)^(-2)*sum((1./model.nRepeats{k}).*...
                                    (model.Dinv{k}.*(model.Ktilde{k} - (model.D{k} - model.m{k}.*model.m{k} + H)).*model.Dinv{k}));
                            end                                                              
                        else
                            error('Information about number of repetitions is not provided')
                        end
                    case 2
                        error('Not implemented yet')                        
                end
            else
                for r =1:model.nlf,
                    for k= 1:model.nout,
                         dLdKuy{r,k} = model.KuuinvKuy{r,k}*model.beta(k) - CKuy{r,k}*model.beta(k) + ...
                            model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);                       
                    end
                end
                for k =1:model.nout,
                    H = zeros(size(model.X{k+model.nlf},1),1);
                    for r =1:model.nlf,
                        temp =  model.Kyu{k,r}*model.AinvKuyDinvy{r};
                        H = H + 2*temp.*model.m{k} - sum(model.Kyu{k,r}.*CKuy{r,k}',2);
                    end
                     dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*1/model.beta(k)...
                        - model.m{k}'*model.m{k} + sum(H)));
                end
                if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)     
                    for k =1:model.nout,
                        DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
                        for r =1:model.nlf,
                            DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
                                model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
                        end
                        dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                    end
                end
            end
            KuuinvKuyQ = cell(model.nlf,model.nout);
            for r =1:model.nlf,
                for k= 1:model.nout,
                    tempoD = model.Dinv{k}(:,ones(1,model.k(r)));
                    KuuinvKuyQ{r,k} = model.KuuinvKuy{r,k}.*tempoD';
                end
            end
            dLdKuu = cell(model.nlf,1);
            for r =1:model.nlf,
                KuuinvKuyQKyuKuuinv = zeros(model.k(r));
                dLdKuu{r} = zeros(model.k(r));
                for k=1: model.nout,
                    KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv ...
                        + KuuinvKuyQ{r,k}*model.KuuinvKuy{r,k}';
                end
                dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
            end
%             dLdKuy = cell(model.nlf,model.nout);
%             for r =1:model.nlf,
%                 for k= 1:model.nout,
%                     if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                         switch model.noiseOpt
%                             case 0
%                                 dLdKuy{r,k} = model.KuuinvKuy{r,k}*model.beta(k)...
%                                     - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
%                             case 1
%                                 if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
%                                     dLdKuy{r,k} = model.KuuinvKuy{r,k}*diag(model.beta(k)*model.nRepeats{k})...
%                                         - CKuy{r,k}*diag(model.beta(k)*model.nRepeats{k}) ...
%                                         + model.AinvKuyDinvy{r}*model.m{k}'*diag(model.beta(k)*model.nRepeats{k});
% 
%                                 else
%                                     error('Information about number of repetitions is not provided')
%                                 end
%                             case 2
%                                 dLdKuy{r,k} = model.KuuinvKuy{r,k}*diag(model.beta{k}) - CKuy{r,k}*diag(model.beta{k})...
%                                     + model.AinvKuyDinvy{r}*model.m{k}'*diag(model.beta{k});
%                         end
%                     else
%                         dLdKuy{r,k} = model.KuuinvKuy{r,k}*model.beta(k)...
%                             - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
%                     end
%                 end
%             end
%             dLdKuu = cell(model.nlf,1);
%             for r =1:model.nlf,
%                 KuuinvKuyQKyuKuuinv = zeros(model.k(r));
%                 dLdKuu{r} = zeros(model.k(r));
%                 for k=1: model.nout,
%                     KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + ...
%                         model.KuuinvKuy{r,k}*model.Dinv{k}*model.KuuinvKuy{r,k}';
%                 end
%                 dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
%             end
%             dLdbeta = zeros(1,model.nout);
%             dLdmu = cell(model.nout,1);
%             for k =1:model.nout,
%                 DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
%                 H = zeros(size(model.X{k+model.nlf},1));
%                 for r =1:model.nlf,
%                     H = H + model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}'+ ...
%                         (model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}')'...
%                         - model.Kyu{k,r}*CKuy{r,k};
%                     DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + ...
%                         model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
%                 end
%                 if ~strcmp(model.kernType,'gg')
%                     if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                         switch model.noiseOpt
%                             case 0
%                                 dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
%                             case 1
%                                 dLdmu{k} = (model.beta(k)*model.nRepeats{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
%                             case 2
%                                 dLdmu{k} = (model.beta{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
%                         end
%                     else
%                         dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
%                     end
%                 end
%                 if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
%                     switch model.noiseOpt
%                         case 0
%                             dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*1/model.beta(k) - ...
%                                 model.m{k}'*model.m{k} + trace(H)));
%                         case 1
%                             dLdbeta(k) = -0.5*model.beta(k)^(-2)*trace(sparseDiag((1./model.nRepeats{k}))*...
%                                 (model.Dinv{k}*(sparseDiag(model.Ktilde{k}) - ...
%                                 (model.D{k} - model.m{k}*model.m{k}' + H))*model.Dinv{k}));
%                         case 2
%                             error('Not implemented yet')
%                     end
%                 else
%                     dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*...
%                         1/model.beta(k) - model.m{k}'*model.m{k} + trace(H))  );
%                 end
%             end
    end
end
