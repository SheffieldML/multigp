function ll = spmultigpLogLikelihood( model)

% SPMULTIGPLOGLIKELIHOOD Compute the log likelihood of a SPARSE MULTIGP.
% FORMAT
% DESC log likelihood in the sparse multigp model
% ARG model : the model structure containing the sparse model
% RETURN ll : log likelihood for the sparse multigp model
%
%
% COPYRIGHT : Mauricio A Alvarez, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2010

% MULTIGP



ll = model.N*log(2*pi);
ll = ll + model.logDetDT; % contribution of ln det D
ll = ll + model.traceDinvyy; % contribution of trace(inv D yy')
ll = ll - model.logDetKuuT; % contribution of ln det Kuu
ll = ll + model.logDetA; % contribution of ln det A
for k =1: model.nlf, % contribution of trace(invD Kyu invA Kuy invD yy')
   ll = ll - sum(model.sqrtAinvKuyDinvy{k}.*model.sqrtAinvKuyDinvy{k});
end
ll = ll + model.KtildeT; % Only meaningful if the approximation is DTCVAR
ll = -0.5*ll;

% %     for k = 1: model.nout,
% %         ll = ll + model.logDetD{k}; % contribution of ln det D
% %         switch model.approx
% %             case {'dtc', 'dtcvar'}
% %                 if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
% %                     switch model.noiseOpt
% %                         case 0
% %                             ll = ll + model.beta(k)*sum(model.m{k}.^2); % contribution of trace(inv D yy')
% %                         case 1
% %                             if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
% %                                 ll = ll + sum(model.beta(k)*model.nRepeats{k}.*(model.m{k}.^2)); % contribution of trace(inv D yy')
% %                             end
% %                         case 2
% %                         case 3
% %                             ll = ll + sum((model.beta{k}).*(model.m{k}.^2)); % contribution of trace(inv D yy')
% %                     end
% %                 else
% %                     ll = ll + model.beta(k)*sum(model.m{k}.^2); % contribution of trace(inv D yy')
% %                 end
% %             case {'fitc','pitc'}
% %                 ll = ll + trace(model.Dinv{k}*model.m{k}*model.m{k}'); % contribution of trace(inv D yy')
% %         end
% %     end
% %     for r = 1:model.nlf, % contribution of ln det Kuu
% %         ll = ll - model.logDetKuu{r};
% %     end
% %     ll = ll + model.logDetA; % contribution of ln det A
% %     for k =1: model.nlf, % contribution of trace(invD Kyu invA Kuy invD yy')
% %         for r = 1:model.nlf,
% %             ll = ll - trace(model.KuyDinvy{k}'*model.Ainv{k,r}*model.KuyDinvy{r});
% %         end
% %     end
% %     if strcmp(model.approx,'dtcvar')
% %         for k = 1:model.nout
% %             if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
% %                 switch model.noiseOpt
% %                     case 0
% %                         ll = ll + model.beta(k)*sum(model.Ktilde{k});
% %                     case 1
% %                         if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
% %                             ll = ll + sum(model.beta(k)*model.nRepeats{k}.*model.Ktilde{k});
% %                         end
% %                     case 2
% %
% %                     case 3
% %                         ll = ll + sum((model.beta{k}).*model.Ktilde{k});
% %                 end
% %             else
% %                 ll = ll + model.beta(k)*sum(model.Ktilde{k});
% %             end
% %         end
% %     end
%     ll = -0.5*ll;


