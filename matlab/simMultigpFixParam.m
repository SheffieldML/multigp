function model = simMultigpFixParam(model, options)

% SIMMULTIGPFIXPARAM Fix the parameters for a multigp model with SIM kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses SIM kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP

if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    % This code fixes latent variances to 1.
    index = paramNameRegularExpressionLookup(model, ['multi ' ...
        '[0-' num2str(options.typeLf(1)) ']+ rbf [0-9]+ variance']);
    count = 0;
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1, 'xtoa');
    end
    % This code fixes the sim-white latent variances to (1/2)^i (i=0, 1, 2, ..., nlf-1).
    for i = 1+options.typeLf(1):model.nlf
        index = paramNameRegularExpressionLookup(model, ...
            ['multi [' num2str(i) '] .* variance']);
        for k=1:length(index);
            count = count + 1;
            model.fix(count).index = index(k);
            model.fix(count).value = expTransform(2^(-i+options.typeLf(1)+1), 'xtoa');
        end
    end
    % This code fixes the variance of the noise in the latent forces to 1e-9.
    for i = 1:model.nlf
        index = paramNameRegularExpressionLookup(model, ...
            ['multi [' num2str(model.nlf+1) '] white [' num2str(i) '] .*']);
        for k=1:length(index);
            count = count + 1;
            model.fix(count).index = index(k);
            model.fix(count).value = expTransform(1e-9, 'xtoa');
        end
    end
else
    % This code fixes masses and latent variances to 1.
    if strcmp(model.kernType, 'sim'),
        index = paramNameRegularExpressionLookup(model, ['multi ' ...
            '[0-' num2str(model.nlf) ']+ rbf [0-9]+ variance']);
    else
        index = paramNameRegularExpressionLookup(model, ['multi ' ...
            '[0-' num2str(model.nlf) ']+ rbfnorm [0-9]+ variance']);
    end
    count = 0;
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1, 'xtoa');
    end
    if ~strcmp(model.approx, 'ftc')
        if options.includeNoise
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
        end
        if isfield(options, 'connect') && ~isempty(options.connect)         
            indexSens = paramNameRegularExpressionLookup(model, '.* sim .* variance');
            temp = options.connect(:);
            indexZeros = find(temp == 0);            
            for k =1:length(indexZeros),
                count = count + 1;
                model.fix(count).index = indexSens(indexZeros(k));
                model.fix(count).value = expTransform(0, 'xtoa');
            end                        
        end
    end
end