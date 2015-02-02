function model = multigpUpdateKernels(model)

% MULTIGPUPDATEKERNELS Updates the kernel representations in the MULTIGP structure.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if isfield(model, 'isSpeedUp') && ~isempty(model.isSpeedUp) && (model.isSpeedUp == 2)
    if strcmp(model.approx, 'ftc')
        K = globalKernCompute(model.kern, model.X(model.nlf+1:end));
        if ~model.includeNoise && isfield(model, 'beta') && ~isempty(model.beta)            
            startVal = 1;
            endVal = 0;
            for k=1:model.nout
                endVal = endVal + size(model.X{model.nlf+k},1);
                K(startVal:endVal,startVal:endVal) = K(startVal:endVal,startVal:endVal) ...
                    + diag((1/model.beta(k))*(ones(size(model.X{model.nlf+k},1),1)));
                startVal = endVal + 1;
            end            
        end
        model.K = K;
%        model.invK = inv(K);
        [model.invK, U, jitter] = pdinv(K);
        if any(jitter>1e-4)
            fprintf('Warning: UpdateKernels added jitter of %2.4f\n', jitter)
        end
        model.alpha = multigpComputeAlpha(model);
        model.logDetK = logdet(model.K, U);
    else
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            [model.Kyy, model.Kyu, model.KuuGamma] = globalKernCompute(model.kern, ...
                model.X(model.nlf+1:end), model.X(1:model.nlf), model.gamma);
        else
            [model.Kyy, model.Kyu, model.Kuu] = globalKernCompute(model.kern, ...
                model.X(model.nlf+1:end), model.X(1:model.nlf));
        end
        model = spmultigpUpdateAD(model);
    end
else
    switch model.approx
        case 'ftc'
            if strcmp(model.kernType, 'lmc')
                K = zeros(model.N);
                for i=1:length(model.kern.comp)
                    if model.isIsotopic
                        K = K + kernCompute(model.kern.comp{i}.comp{model.nlf+1}, ....
                            model.X{model.nlf+1});
                    else
                        K = K + kernCompute(model.kern.comp{i}.comp{model.nlf+1}, ....
                            model.X(model.nlf+1:end));
                    end
                end
            else
                K = real(kernCompute(model.kern, model.X));
                % This is a hack to allow the inclusion of the latent structure
                % inside the kernel structure
                if sum(cellfun('length',model.X)) == model.N,
                    base = 0;
                else
                    base = model.nlf;
                end
                K = K(base+1:end, base+1:end);
            end
            
            if ~model.includeNoise && isfield(model, 'beta') && ~isempty(model.beta)
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt) && model.noiseOpt
                    startVal = 1;
                    endVal = 0;
                    for k=1:model.nout
                        endVal = endVal + size(model.X{model.nlf+k},1);
                        K(startVal:endVal,startVal:endVal) = K(startVal:endVal,startVal:endVal) ...
                            + diag((1/model.beta(k))*(1./model.nRepeats{k}));
                        startVal = endVal + 1;
                    end
                end
            end
            model.K = K;
            %/~
            %     [L, jitter] = blockChol(model);
            %     invL = L\eye(size(K,1));
            %     model.invK = invL'*invL;
            %~/
            model.invK = inv(K);
            [model.invK, U, jitter] = pdinv(K);
            model.logDetK = logdet(model.K, U);
            %
            %L = chol(K);
            %Linv = L\eye(size(L));
            %model.invK  = Linv*Linv';
            %model.logDetK = 2*sum(log(diag(L)));
            %
            model.alpha = multigpComputeAlpha(model);
            %/~
            %model.logDetK = 2*sum(log(diag(L)));
            %~/                        
        case {'dtc','fitc','pitc', 'dtcvar'}
            model = spmultigpUpdateKernels(model);
    end
end

