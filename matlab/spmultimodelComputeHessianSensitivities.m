function [Sens, gS, gS2] = spmultimodelComputeHessianSensitivities(model, tf, mode)

% SPMULTIMODELCOMPUTEHESSIANSENSITIVITIES Laplace appproximation for LFM.
% Computes the Laplace approximation for the sensitivity parameters of a
% latent force model with kernel 'SIM'. This was particularly made for the
% cell cycle Yeast example due to Spellman et. al. (1998).
% FORMAT
% DESC returns sensitivities of the model, the gradient and the second
% derivative of the sensitivities associated to a particular latent force.
% RETURN Sens : value of the sensitivity parameters
% RETURN gS   : first derivative of the sensitivity parameters
% RETURN Sens : second derivative of the sensitivity parameters
% ARG model : the 'spmultimodel' structure.
% ARG tf : the number of the transcription factor whose associated
% sensitivities we are interested in.
% ARG mode : the mode is the method employed to compute the gradients. It
% can be 'numeric' in which case the derivatives are computed numerically
% or it can be 'computeTest' in which case the derivatives are computed
% analitically. An additional option is 'computeTestDecay' which computes
% the gradients of the parameter defined as the sensitivity over the decay
% in the SIM kernel.
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP


Sens = zeros(1, model.numOutGivenLat(tf));
gS = zeros(1, model.numOutGivenLat(tf));
gS2 = zeros(model.numOutGivenLat(tf));

if nargin < 3
    % Using the parameters of the sparse in the full model
    mode = 'computeTest';
end


for i=1:model.numOutGivenLat(tf)
    Sens(i) = model.kern.comp{tf}.comp{i+1}.variance;
end

switch mode
    case 'numeric'
        % Compute numerically the second derivative
        subIndex = paramNameRegularExpressionLookup(model, ['component ' num2str(tf) ' .* variance']);
        subIndex = subIndex(2:end);
        epsilon = 1e-4;
        paramsBase = modelExtractParam(model);
        %params = paramsBase(subIndex);
        %         w = paramsBase;
        %         func = 'modelObjective';
        %         nparams = length(subIndex);
        %         stepi = zeros(1, length(paramsBase));
        %         stepj = zeros(1, length(paramsBase));
        %         f = feval('linef', epsilon, func, w, stepi, model);
        %         for i = 1:nparams
        %             stepi(subIndex(i)) = 1.0;
        %             fi = feval('linef', epsilon, func, w, stepi, model);
        %             fj = feval('linef', epsilon, func, w, stepi, model);
        %             fij = feval('linef', epsilon, func, w, stepi+stepi, model);
        %             gS2(i, i) = (fij- fi - fj + f)/(epsilon^2);
        %             for j=1:i-1,
        %                 stepj(subIndex(j)) = 1.0;
        %                 fj = feval('linef', epsilon, func, w, stepj, model);
        %                 fij = feval('linef', epsilon, func, w, stepi+stepj, model);
        %                 gS2(i, j) = (fij- fi - fj + f)/(epsilon^2);
        %                 stepj(subIndex(j)) = 0;
        %             end
        %             stepi(subIndex(i)) = 0;
        %         end
        %         for i = 1:nparams
        %             for j=1:i-1
        %                 gS2(j,i) = gS2(i,j);
        %             end
        %         end
        params = paramsBase(subIndex);
        origParams = params;
        for i = 1:length(params);
            params = origParams;
            params(i) = origParams(i) + epsilon;
            paramsBase(subIndex) = params;
            model = modelExpandParam(model, paramsBase);
            Lplus = modelLogLikelihood(model);
            params(i) = origParams(i) - epsilon;
            paramsBase(subIndex) = params;
            model = modelExpandParam(model, paramsBase);
            Lminus = modelLogLikelihood(model);
            gS(i) = 0.5*(Lplus - Lminus)/epsilon;
            paramsBase(subIndex) = origParams;
            model = modelExpandParam(model, paramsBase);
            L = modelLogLikelihood(model);
            gS2(i) = (Lplus - 2*L + Lminus)/(epsilon^2);
        end
    case 'computeTest'

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Kffinv for all the outputs for a test of the derivatives in small datasets
        %
        %         C = model.connect*model.latConnectInv;
        %         Kffinv = cell(model.nout);
        %         for n = 1:model.nout,
        %             whichLatent = model.indLatGivenOut{n};
        %             KfuAinvKuf = zeros(model.sizeX(n));
        %             whichColumnAinv = find(C(n, :)~=0);
        %             KfuAinv = cell(length(whichColumnAinv),1);
        %             for r=1:length(whichColumnAinv)
        %                 KfuAinv{r} = zeros(model.sizeX(n), model.k);
        %                 whichLatent = model.indLatGivenOut{n};
        %                 for s = 1:model.numLatGivenOut(n),
        %                     KfuAinv{r} = KfuAinv{r} + ...
        %                         model.Kyu{n}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
        %                 end
        %             end
        %             for r =1:model.numLatGivenOut(n)
        %                 KfuAinvKuf = KfuAinvKuf + ...
        %                     KfuAinv{whichColumnAinv ==whichLatent(r)}*(model.Kyu{n}{r})';
        %             end
        %             Kffinv{n,n} = model.Dinv{n} ...
        %                 - model.Dinv{n}*KfuAinvKuf*model.Dinv{n};
        %             for m = 1:n-1
        %                 whichLatent2 = model.indLatGivenOut{m};
        %                 KfuAinvKuf = zeros(model.sizeX(n));
        %                 for r =1:model.numLatGivenOut(m)
        %                     KfuAinvKuf = KfuAinvKuf + ...
        %                         KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{m}{r})';
        %                 end
        %                 Kffinv{n,m} = - model.Dinv{n}*KfuAinvKuf*model.Dinv{m};
        %                 Kffinv{m,n} = (Kffinv{n,m})';
        %             end
        %         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Just to compute the bit of the Kffinv that it's necessary.

        C = model.connect*model.latConnectInv;
        whichOutputsGlobal = model.indOutGivenLat{tf};
        Kffinv = cell(model.numOutGivenLat(tf));
        for n = 1:model.numOutGivenLat(tf)
            whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
            KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
            whichColumnAinv = find(C(whichOutputsGlobal(n), :)~=0);
            KfuAinv = cell(length(whichColumnAinv),1);
            for r=1:length(whichColumnAinv)
                KfuAinv{r} = zeros(model.sizeX(whichOutputsGlobal(n)), model.k);
                whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
                for s = 1:model.numLatGivenOut(whichOutputsGlobal(n)),
                    KfuAinv{r} = KfuAinv{r} + ...
                        model.Kyu{whichOutputsGlobal(n)}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
                end
            end
            for r =1:model.numLatGivenOut(whichOutputsGlobal(n))
                KfuAinvKuf = KfuAinvKuf + ...
                    KfuAinv{whichColumnAinv ==whichLatent(r)}*(model.Kyu{whichOutputsGlobal(n)}{r})';
            end
            Kffinv{n,n} = model.Dinv{whichOutputsGlobal(n)} ...
                - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(n)};
            for m = 1:n-1
                whichLatent2 = model.indLatGivenOut{whichOutputsGlobal(m)};
                KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
                for r =1:model.numLatGivenOut(whichOutputsGlobal(m))
                    KfuAinvKuf = KfuAinvKuf + ...
                        KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{whichOutputsGlobal(m)}{r})';
                end
                Kffinv{n,m} = - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(m)};
                Kffinv{m,n} = (Kffinv{n,m})';
            end
        end
        KffinvSquare = Kffinv;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Just to compute the bit of Kffinv*m that it's necessary
        % First the lower part

        C2 = pdinv(model.connect*model.connect' + 1e3*eye(size(model.connect,1)));
        C = model.connect*model.latConnectInv;
        whichOutputsGlobal = model.indOutGivenLat{tf};
        Kffinv = cell(model.numOutGivenLat(tf));
        for n = 1:model.numOutGivenLat(tf)
            whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
            KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
            whichColumnAinv = find(C(whichOutputsGlobal(n), :)~=0);
            KfuAinv = cell(length(whichColumnAinv),1);
            for r=1:length(whichColumnAinv)
                KfuAinv{r} = zeros(model.sizeX(whichOutputsGlobal(n)), model.k);
                whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
                for s = 1:model.numLatGivenOut(whichOutputsGlobal(n)),
                    KfuAinv{r} = KfuAinv{r} + ...
                        model.Kyu{whichOutputsGlobal(n)}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
                end
            end
            for r =1:model.numLatGivenOut(whichOutputsGlobal(n))
                KfuAinvKuf = KfuAinvKuf + ...
                    KfuAinv{whichColumnAinv ==whichLatent(r)}*(model.Kyu{whichOutputsGlobal(n)}{r})';
            end
            Kffinv{n,whichOutputsGlobal(n)} = model.Dinv{whichOutputsGlobal(n)} ...
                - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(n)};
            for m = 1:whichOutputsGlobal(n)-1
                if C2(whichOutputsGlobal(n),m)~=0
                    whichLatent2 = model.indLatGivenOut{m};
                    KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
                    for r =1:model.numLatGivenOut(m)
                        KfuAinvKuf = KfuAinvKuf + ...
                            KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{m}{r})';
                    end
                    Kffinv{n,m} = - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{m};
                else
                    Kffinv{n,m} = [];
                end
            end
        end
        % The the upper part
        for n = 1:model.numOutGivenLat(tf)
            for m =whichOutputsGlobal(n)+1:model.nout
                if C2(m, whichOutputsGlobal(n))~=0
                    KfuAinvKuf = zeros(model.sizeX(m));
                    whichColumnAinv = find(C(m, :)~=0);
                    KfuAinv = cell(length(whichColumnAinv),1);
                    for r=1:length(whichColumnAinv)
                        KfuAinv{r} = zeros(model.sizeX(m), model.k);
                        whichLatent = model.indLatGivenOut{m};
                        for s = 1:model.numLatGivenOut(m),
                            KfuAinv{r} = KfuAinv{r} + ...
                                model.Kyu{m}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
                        end
                    end
                    whichLatent2 = model.indLatGivenOut{whichOutputsGlobal(n)};
                    for r =1:model.numLatGivenOut(whichOutputsGlobal(n))
                        KfuAinvKuf = KfuAinvKuf + ...
                            KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{whichOutputsGlobal(n)}{r})';
                    end
                    Kffinv{n,m} = (- model.Dinv{m}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(n)})';
                else
                    Kffinv{n,m} = [];
                end
            end
        end
        KffinvRect = Kffinv;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply for y
        Kffinvy = cell(model.numOutGivenLat(tf),1);
        for n = 1:model.numOutGivenLat(tf)
            Kffinvy{n} = zeros(model.sizeX(whichOutputsGlobal(n)),1);
            for m=1:model.nout
                if C2(whichOutputsGlobal(n), m)~=0
                    Kffinvy{n} = Kffinvy{n} + KffinvRect{n,m}*model.m{m};
                end
            end
        end
        Kffinv = cell2mat(KffinvSquare);
        KffinvM  = cell2mat(Kffinvy);
        KffinvM2 = KffinvM*KffinvM';
        Qfuf = cell(model.numOutGivenLat(tf));
        Kfftilde = cell(model.numOutGivenLat(tf));
        whichOutputs = model.indOutGivenLat{tf};
        for i=1:model.numOutGivenLat(tf)
            [k, sk1] = simXrbfKernCompute(model.kern.comp{tf}.comp{i+1},...
                model.kern.comp{tf}.comp{1}, model.outX{whichOutputs(i)},...
                model.latX{tf});
            if ~model.kern.comp{tf}.comp{i+1}.isNormalised
                sigma = sqrt(2/model.kern.comp{tf}.comp{i+1}.inverseWidth);
                sk1 = sigma*sk1;
            end
            [k, skf1] = simKernCompute(model.kern.comp{tf}.comp{i+1},...
                model.outX{whichOutputs(i)});
            if ~model.kern.comp{tf}.comp{i+1}.isNormalised
                sigma = sqrt(2/model.kern.comp{tf}.comp{i+1}.inverseWidth);
                skf1 = sqrt(pi)*sigma*skf1;
            end
            switch model.approx
                case {'dtc', 'dtcvar'}
                    Qfuf{i,i} = sk1*model.Kuuinv{tf}*sk1';
                case {'pitc'}
                    Qfuf{i,i} = skf1;
            end
            Kfftilde{i,i} = skf1;
            for j = 1:i-1
                [k, sk2] = simXrbfKernCompute(model.kern.comp{tf}.comp{j+1},...
                    model.kern.comp{tf}.comp{1}, model.outX{whichOutputs(j)},...
                    model.latX{tf});
                if ~model.kern.comp{tf}.comp{j+1}.isNormalised
                    sigma = sqrt(2/model.kern.comp{tf}.comp{j+1}.inverseWidth);
                    sk2 = sigma*sk2;
                end
                [k, skf2] = simXsimKernCompute(model.kern.comp{tf}.comp{i+1},...
                    model.kern.comp{tf}.comp{j+1}, model.outX{whichOutputs(i)},...
                    model.outX{whichOutputs(j)});
                if ~model.kern.comp{tf}.comp{j+1}.isNormalised
                    sigma = sqrt(2/model.kern.comp{tf}.comp{j+1}.inverseWidth);
                    skf2 = sqrt(pi)*sigma*skf2;
                end
                Qfuf{i,j} = sk1*model.Kuuinv{tf}*sk2';
                Qfuf{j,i} = Qfuf{i,j}';
                Kfftilde{i,j} = skf2;
                Kfftilde{j,i} = skf2';
            end
        end
        for i=1:model.numOutGivenLat(tf),
            for j=1:model.numOutGivenLat(tf),
                if isempty(Kfftilde{i,j})
                    Qfuf{i, j} = zeros(model.sizeX(i),model.sizeX(j));
                    Kfftilde{i, j} = zeros(model.sizeX(i),model.sizeX(j));
                end
            end
        end
        % For the S matrices according to the needed parameter
        QfufF = cell2mat(Qfuf);
        KfftildeF = cell2mat(Kfftilde);
        for i = 1:model.numOutGivenLat(tf)
            S = zeros(model.numOutGivenLat(tf));
            S2 = zeros(model.numOutGivenLat(tf));
            for j =1:model.numOutGivenLat(tf)
                for k=1:model.numOutGivenLat(tf)
                    if j==i || k==i
                        if j==i && k==i
                            S(j,j) =  2*model.kern.comp{tf}.comp{1+j}.variance;
                            S2(j,j) = 2;
                        elseif j==i && k~=i
                            S(j,k) = model.kern.comp{tf}.comp{1+k}.variance;
                        elseif j~=i && k==i
                            S(j,k) = model.kern.comp{tf}.comp{1+j}.variance;
                        end
                    end
                end
            end
            bigS = kron(S, ones(model.sizeX(1)));
            bigS2 = kron(S2, ones(model.sizeX(1)));
            dKffdS = bigS.*QfufF;
            %%% To get the first derivative
            gS(i) = -0.5*trace((Kffinv - KffinvM2)'*dKffdS);
            %%% To get second derivative
            partialdKff2 = Kffinv*dKffdS*Kffinv - Kffinv*dKffdS*KffinvM2 - ...
                KffinvM2*dKffdS*Kffinv;
            gS2(i,i) = 0.5*(partialdKff2(:)'*dKffdS(:))- ...
                0.5*trace((Kffinv- KffinvM2)'*(bigS2.*QfufF));
            if strcmp(model.approx, 'dtcvar')
                beta = (1./model.beta(model.indOutGivenLat{tf}))';
                beta = kron(beta, ones(size(model.outX{1},1),1));
                gS(i) = gS(i) - 0.5*trace(diag(1./beta)*(bigS.*(KfftildeF - QfufF)));
                gS2(i,i) = gS2(i,i) - 0.5*trace(diag(1./beta)*(bigS2.*(KfftildeF - QfufF)));
            end
        end

    case 'computeTestDecay'

        % Just to compute the bit of the Kffinv that it's necessary.

        C = model.connect*model.latConnectInv;
        whichOutputsGlobal = model.indOutGivenLat{tf};
        Kffinv = cell(model.numOutGivenLat(tf));
        for n = 1:model.numOutGivenLat(tf)
            whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
            KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
            whichColumnAinv = find(C(whichOutputsGlobal(n), :)~=0);
            KfuAinv = cell(length(whichColumnAinv),1);
            for r=1:length(whichColumnAinv)
                KfuAinv{r} = zeros(model.sizeX(whichOutputsGlobal(n)), model.k);
                whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
                for s = 1:model.numLatGivenOut(whichOutputsGlobal(n)),
                    KfuAinv{r} = KfuAinv{r} + ...
                        model.Kyu{whichOutputsGlobal(n)}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
                end
            end
            for r =1:model.numLatGivenOut(whichOutputsGlobal(n))
                KfuAinvKuf = KfuAinvKuf + ...
                    KfuAinv{whichColumnAinv ==whichLatent(r)}*(model.Kyu{whichOutputsGlobal(n)}{r})';
            end
            Kffinv{n,n} = model.Dinv{whichOutputsGlobal(n)} ...
                - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(n)};
            for m = 1:n-1
                whichLatent2 = model.indLatGivenOut{whichOutputsGlobal(m)};
                KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
                for r =1:model.numLatGivenOut(whichOutputsGlobal(m))
                    KfuAinvKuf = KfuAinvKuf + ...
                        KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{whichOutputsGlobal(m)}{r})';
                end
                Kffinv{n,m} = - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(m)};
                Kffinv{m,n} = (Kffinv{n,m})';
            end
        end
        KffinvSquare = Kffinv;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Just to compute the bit of Kffinv*m that it's necessary
        % First the lower part

        C2 = pdinv(model.connect*model.connect' + 1e3*eye(size(model.connect,1)));
        C = model.connect*model.latConnectInv;
        whichOutputsGlobal = model.indOutGivenLat{tf};
        Kffinv = cell(model.numOutGivenLat(tf));
        for n = 1:model.numOutGivenLat(tf)
            whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
            KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
            whichColumnAinv = find(C(whichOutputsGlobal(n), :)~=0);
            KfuAinv = cell(length(whichColumnAinv),1);
            for r=1:length(whichColumnAinv)
                KfuAinv{r} = zeros(model.sizeX(whichOutputsGlobal(n)), model.k);
                whichLatent = model.indLatGivenOut{whichOutputsGlobal(n)};
                for s = 1:model.numLatGivenOut(whichOutputsGlobal(n)),
                    KfuAinv{r} = KfuAinv{r} + ...
                        model.Kyu{whichOutputsGlobal(n)}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
                end
            end
            for r =1:model.numLatGivenOut(whichOutputsGlobal(n))
                KfuAinvKuf = KfuAinvKuf + ...
                    KfuAinv{whichColumnAinv ==whichLatent(r)}*(model.Kyu{whichOutputsGlobal(n)}{r})';
            end
            Kffinv{n,whichOutputsGlobal(n)} = model.Dinv{whichOutputsGlobal(n)} ...
                - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(n)};
            for m = 1:whichOutputsGlobal(n)-1
                if C2(whichOutputsGlobal(n),m)~=0
                    whichLatent2 = model.indLatGivenOut{m};
                    KfuAinvKuf = zeros(model.sizeX(whichOutputsGlobal(n)));
                    for r =1:model.numLatGivenOut(m)
                        KfuAinvKuf = KfuAinvKuf + ...
                            KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{m}{r})';
                    end
                    Kffinv{n,m} = - model.Dinv{whichOutputsGlobal(n)}*KfuAinvKuf*model.Dinv{m};
                else
                    Kffinv{n,m} = [];
                end
            end
        end
        % The the upper part
        for n = 1:model.numOutGivenLat(tf)
            for m =whichOutputsGlobal(n)+1:model.nout
                if C2(m, whichOutputsGlobal(n))~=0
                    KfuAinvKuf = zeros(model.sizeX(m));
                    whichColumnAinv = find(C(m, :)~=0);
                    KfuAinv = cell(length(whichColumnAinv),1);
                    for r=1:length(whichColumnAinv)
                        KfuAinv{r} = zeros(model.sizeX(m), model.k);
                        whichLatent = model.indLatGivenOut{m};
                        for s = 1:model.numLatGivenOut(m),
                            KfuAinv{r} = KfuAinv{r} + ...
                                model.Kyu{m}{s}*model.Ainv{whichLatent(s),whichColumnAinv(r)};
                        end
                    end
                    whichLatent2 = model.indLatGivenOut{whichOutputsGlobal(n)};
                    for r =1:model.numLatGivenOut(whichOutputsGlobal(n))
                        KfuAinvKuf = KfuAinvKuf + ...
                            KfuAinv{whichColumnAinv ==whichLatent2(r)}*(model.Kyu{whichOutputsGlobal(n)}{r})';
                    end
                    Kffinv{n,m} = (- model.Dinv{m}*KfuAinvKuf*model.Dinv{whichOutputsGlobal(n)})';
                else
                    Kffinv{n,m} = [];
                end
            end
        end
        KffinvRect = Kffinv;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply for y
        Kffinvy = cell(model.numOutGivenLat(tf),1);
        for n = 1:model.numOutGivenLat(tf)
            Kffinvy{n} = zeros(model.sizeX(whichOutputsGlobal(n)),1);
            for m=1:model.nout
                if C2(whichOutputsGlobal(n), m)~=0
                    Kffinvy{n} = Kffinvy{n} + KffinvRect{n,m}*model.m{m};
                end
            end
        end
        Kffinv = cell2mat(KffinvSquare);
        KffinvM  = cell2mat(Kffinvy);
        KffinvM2 = KffinvM*KffinvM';
        Qfuf = cell(model.numOutGivenLat(tf));
        Kfftilde = cell(model.numOutGivenLat(tf));
        whichOutputs = model.indOutGivenLat{tf};
        for i=1:model.numOutGivenLat(tf)
            [k, sk1] = simXrbfKernCompute(model.kern.comp{tf}.comp{i+1},...
                model.kern.comp{tf}.comp{1}, model.outX{whichOutputs(i)},...
                model.latX{tf});
            if ~model.kern.comp{tf}.comp{i+1}.isNormalised
                sigma = sqrt(2/model.kern.comp{tf}.comp{i+1}.inverseWidth);
                sk1 = sigma*sk1;
            end
            [k, skf1] = simKernCompute(model.kern.comp{tf}.comp{i+1},...
                model.outX{whichOutputs(i)});
            if ~model.kern.comp{tf}.comp{i+1}.isNormalised
                sigma = sqrt(2/model.kern.comp{tf}.comp{i+1}.inverseWidth);
                skf1 = sqrt(pi)*sigma*skf1;
            end
            switch model.approx
                case {'dtc', 'dtcvar'}
                    Qfuf{i,i} = sk1*model.Kuuinv{tf}*sk1';
                case {'pitc'}
                    Qfuf{i,i} = skf1;
            end
            Kfftilde{i,i} = skf1;
            for j = 1:i-1
                [k, sk2] = simXrbfKernCompute(model.kern.comp{tf}.comp{j+1},...
                    model.kern.comp{tf}.comp{1}, model.outX{whichOutputs(j)},...
                    model.latX{tf});
                if ~model.kern.comp{tf}.comp{j+1}.isNormalised
                    sigma = sqrt(2/model.kern.comp{tf}.comp{j+1}.inverseWidth);
                    sk2 = sigma*sk2;
                end
                [k, skf2] = simXsimKernCompute(model.kern.comp{tf}.comp{i+1},...
                    model.kern.comp{tf}.comp{j+1}, model.outX{whichOutputs(i)},...
                    model.outX{whichOutputs(j)});
                if ~model.kern.comp{tf}.comp{j+1}.isNormalised
                    sigma = sqrt(2/model.kern.comp{tf}.comp{j+1}.inverseWidth);
                    skf2 = sqrt(pi)*sigma*skf2;
                end
                Qfuf{i,j} = sk1*model.Kuuinv{tf}*sk2';
                Qfuf{j,i} = Qfuf{i,j}';
                Kfftilde{i,j} = skf2;
                Kfftilde{j,i} = skf2';
            end
        end
        for i=1:model.numOutGivenLat(tf),
            for j=1:model.numOutGivenLat(tf),
                if isempty(Kfftilde{i,j})
                    Qfuf{i, j} = zeros(model.sizeX(i),model.sizeX(j));
                    Kfftilde{i, j} = zeros(model.sizeX(i),model.sizeX(j));
                end
            end
        end
        % For the S matrices according to the needed parameter
        QfufF = cell2mat(Qfuf);
        KfftildeF = cell2mat(Kfftilde);
        for i = 1:model.numOutGivenLat(tf)
            S = zeros(model.numOutGivenLat(tf));
            S2 = zeros(model.numOutGivenLat(tf));
            for j =1:model.numOutGivenLat(tf)
                for k=1:model.numOutGivenLat(tf)
                    if j==i || k==i
                        if j==i && k==i
                            S(j,j) =  2*model.kern.comp{tf}.comp{1+j}.variance;
                            S2(j,j) = 2;
                        elseif j==i && k~=i
                            S(j,k) = model.kern.comp{tf}.comp{1+k}.variance;
                        elseif j~=i && k==i
                            S(j,k) = model.kern.comp{tf}.comp{1+j}.variance;
                        end
                    end
                end
            end
            S = S*model.kern.comp{tf}.comp{1+j}.decay; % Here is where we put the decay
            S2 = S2*(model.kern.comp{tf}.comp{1+j}.decay^2);% Here is where we put the decay
            bigS = kron(S, ones(model.sizeX(1)));
            bigS2 = kron(S2, ones(model.sizeX(1)));
            dKffdS = bigS.*QfufF;
            %%% To get the first derivative
            gS(i) = -0.5*trace((Kffinv - KffinvM2)'*dKffdS);
            %%% To get second derivative
            partialdKff2 = Kffinv*dKffdS*Kffinv - Kffinv*dKffdS*KffinvM2 - ...
                KffinvM2*dKffdS*Kffinv;
            gS2(i,i) = 0.5*(partialdKff2(:)'*dKffdS(:))- ...
                0.5*trace((Kffinv- KffinvM2)'*(bigS2.*QfufF));
            if strcmp(model.approx, 'dtcvar')
                beta = (1./model.beta(model.indOutGivenLat{tf}))';
                beta = kron(beta, ones(size(model.outX{1},1),1));
                gS(i) = gS(i) - 0.5*trace(diag(1./beta)*(bigS.*(KfftildeF - QfufF)));
                gS2(i,i) = gS2(i,i) - 0.5*trace(diag(1./beta)*(bigS2.*(KfftildeF - QfufF)));
            end
        end
end

