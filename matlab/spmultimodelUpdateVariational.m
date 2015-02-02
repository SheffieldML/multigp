function model = spmultimodelUpdateVariational(model)

% SPMULTIMODELUPDATEVARIATIONAL Updates variational distributions
% FORMAT
% DESC
% Updates the variational distributions over sensitivities and latent
% forces.
% ARG model : model containing the variational distributions to be updated.
% RETURN  model : model with the updated variational distributions.
%
% COPYRIGHT : Mauricio Alvarez, 2010

% MULTIGP

% Update q(u) and compute quantities for q(S)

GammaInv = zeros(model.nlf, model.nlf, model.nout);
aBold = zeros(model.nlf, sum(model.sizeX));
for i=1:model.nlf
    if isfield(model, 'gamma') && ~isempty(model.gamma)
        model.qu.mean{i} = model.KuuGamma{i}*model.AinvKuyHatDinvy{i};
        model.qu.Sigma{i,i} = model.KuuGamma{i}*model.Ainv{i,i}*model.KuuGamma{i};
    else
        model.qu.mean{i} = model.Kuu{i}*model.AinvKuyHatDinvy{i};
        model.qu.Sigma{i,i} = model.Kuu{i}*model.Ainv{i,i}*model.Kuu{i};
    end
    exu2 = model.qu.Sigma{i,i} + model.qu.mean{i}*model.qu.mean{i}';
    Aqq = model.Kuuinv{i}*exu2*model.Kuuinv{i};
    Kuy = cell2mat(model.Kyu(:,i))';
    AqqKuy = Aqq*Kuy;
    GammaS = sum(Kuy.*AqqKuy,1);
    vecGammaS = GammaS*model.template;
%     GammaInv(i,i,:) = model.gammas(((i-1)*model.nout)+1:(i*model.nout)) + ...
%         model.Ed(:,i)' + vecGammaS;
    GammaInv(i,i,:) = model.gammas(((i-1)*model.nout)+1:(i*model.nout)) + ...
        model.beta.*(model.Ed(:,i)' + vecGammaS);
    KuuinvKuy = cell2mat(model.KuuinvKuy(i, :));
    aBold(i,:) = model.qu.mean{i}'*KuuinvKuy;
    for j=1:i-1
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            model.qu.Sigma{i,j} = model.KuuGamma{i}*model.Ainv{i,j}*model.KuuGamma{j};
        else
            model.qu.Sigma{i,j} = model.Kuu{i}*model.Ainv{i,j}*model.Kuu{j};
        end
        model.qu.Sigma{j,i} = model.qu.Sigma{i,j}';
        exu2 = model.qu.Sigma{i,j} + model.qu.mean{i}*model.qu.mean{j}';
        Aqq = model.Kuuinv{i}*exu2*model.Kuuinv{j};
        Kuy2 = cell2mat(model.Kyu(:,j))';
        AqqKuy2 = Aqq*Kuy2;
        GammaS = sum(Kuy.*AqqKuy2,1);
        vecGammaS = GammaS*model.template;
%         GammaInv(i,j,:) = vecGammaS;
%         GammaInv(j,i,:) = vecGammaS;
        GammaInv(i,j,:) = model.beta.*vecGammaS;
        GammaInv(j,i,:) = GammaInv(i,j,:);
    end
end

% Update q(S)
startVal = 1;
endVal = 0;
for i=1:model.nout
    endVal = endVal + model.sizeX(i);
    [model.qs.Sigma(:,:,i), void, jitter] = pdinv(GammaInv(:,:,i));
    yd = repmat(model.beta(i)*model.m{i}', model.nlf, 1);
    aBold2 = yd.*aBold(:,startVal:endVal);
%     model.qs.mean(i,:) = sum(aBold2,2)';
    model.qs.mean(i,:) = (model.qs.Sigma(:,:,i)*sum(aBold2,2))';
    startVal = endVal + 1;
end


