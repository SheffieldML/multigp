function model = ggMultigpFixParam(model, options)

% GGMULTIGPFIXPARAM Fix the parameters for a multigp model with GG kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses GG kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP


if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    index = paramNameRegularExpressionLookup(model, 'multi .* variance latent');
    count = 0;
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1, 'xtoa');
    end
    if (model.nlf - options.typeLf(1)) > 1
        base = 0.5;
    else
        base = 1;
    end
    for i = 1+options.typeLf(1):model.nlf
        index = paramNameRegularExpressionLookup(model, ...
            ['multi [' num2str(i) '] .* variance']);
        for k=1:length(index);
            count = count + 1;
            model.fix(count).index = index(k);
            model.fix(count).value = expTransform(base*5^(i-(1+options.typeLf(1))), 'xtoa');
        end
    end
else
    index = paramNameRegularExpressionLookup(model, 'multi .* variance latent');
    count = 0;
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1, 'xtoa');
    end    
end

if ~strcmp(model.approx, 'ftc')
    % In the approximations fitc, pitc and dtc, this function is
    % accomplished by the parameter beta
    % If there is noise then is at the end and it's not necessary to look
    % for them using paramNameRegularExpressionLookUp
    nParamsKern = 0;
    nToLook = options.nlf + options.includeInd;
    for k =1:nToLook,
        nParamsKern = nParamsKern + model.kern.comp{k}.nParams;        
    end
    index = (nParamsKern + options.nlf + 1):...
        (nParamsKern + options.nlf +model.nout);     
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1e-9, 'xtoa');
    end    
%     for k = 1:model.nout
%         index = paramNameRegularExpressionLookup(model, ...
%             ['multi ' num2str(length(model.kern.comp)) ' white ' ...
%             num2str(model.nlf+k) ' .*']);
%         count = count + 1;
%         model.fix(count).index = index;
%         model.fix(count).value = expTransform(1e-9, 'xtoa');
%     end
end