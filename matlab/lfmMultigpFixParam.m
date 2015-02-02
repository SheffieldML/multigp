function model = lfmMultigpFixParam(model, options)

% LFMMULTIGPFIXPARAM Fix the parameters for a multigp model with LFM kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses LFM kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP

% This code fixes masses and latent variances to 1.
index = paramNameRegularExpressionLookup(model, ['multi ' ...
    '[0-9]+ rbf [0-9]+ variance']);
count = 0;
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = expTransform(1, 'xtoa');
end
% Fix the masses of the lfm kernels.
index = paramNameRegularExpressionLookup(model, ['multi ' ...
    '[0-9]+ lfm [0-9]+ mass']);
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = expTransform(1, 'xtoa');
end

if ~strcmp(model.approx, 'ftc') && options.includeNoise
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