function ll = spmultimodelLogLikelihood( model)

% SPMULTIMODELOGLIKELIHOOD Compute the log likelihood of a SPARSE MULTIMODEL.
% FORMAT
% DESC log likelihood in the sparse multigp model
% ARG model : the model structure containing the sparse model
% RETURN ll : log likelihood for the sparse multigp model
%
% COPYRIGHT : Mauricio A Alvarez, 2009
%
% MODIFICATIONS : Mauricio A Alvarez, 2010

% MULTIGP

ll = sum(model.sizeX)*log(2*pi);

for k = 1: model.nout,
    ll = ll + model.logDetD{k}; % contribution of ln det D
    switch model.approx
        case {'dtc', 'dtcvar'}
            ll = ll + model.beta(k)*sum(model.m{k}.^2); % contribution of trace(inv D yy')
        case {'fitc','pitc'}
            ll = ll + trace(model.Dinv{k}*model.m{k}*model.m{k}'); % contribution of trace(inv D yy')
    end
end
for r = 1:model.nlf, % contribution of ln det Kuu
     ll = ll - model.logDetKuu{r};
end

ll = ll + model.logDetA; % contribution of ln det A


if model.varS
    for k =1: model.nlf, % contribution of trace(invD Kyu invA Kuy invD yy')
        for r = 1:model.nlf,
            ll = ll - model.KuyHatDinvy{k}'*model.Ainv{k,r}*model.KuyHatDinvy{r};
        end
    end
    ll = ll + model.sumCdq; % contribution of c_{d,q}
    ll = ll + model.nout*model.nlf*log(2*pi);    
    ll = ll - sum(log(model.gammas));% contribution of sum_d sum_q log gamma_{d,q} 
    ll = ll + model.sumGdqExS2; % contribution of sum_d sum_q gamma_{d,q}*E[S^2]
    ll = ll - model.entropyS; % contribution of H(S)
else
    for k =1: model.nlf, % contribution of trace(invD Kyu invA Kuy invD yy')
        for r = 1:model.nlf,
            if (model.latConnectInv(k,r)~=0)
                ll = ll - trace(model.KuyDinvy{k}'*model.Ainv{k,r}*model.KuyDinvy{r});
            end
        end
    end
    if strcmp(model.approx,'dtcvar')
        for k = 1:model.nout
            ll = ll + model.beta(k)*sum(model.Ktilde{k});
        end
    end
end

ll = -0.5*ll;
