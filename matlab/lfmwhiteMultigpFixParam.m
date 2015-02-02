function model = lfmwhiteMultigpFixParam(model, options)

% LFMWHITEMULTIGPFIXPARAM Fix parameters for a multigp with LFMWHITE kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses LFMWHITE kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP

% This code fixes the latent variances to 1.
index = paramNameRegularExpressionLookup(model, ...
    ['multi [1-' num2str(model.nlf) '] .* variance']);
count = 0;
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = expTransform(1, 'xtoa');
end
% Fix the masses of the lfm kernels.
index = paramNameRegularExpressionLookup(model, ['multi ' ...
    '[0-9]+ lfmwhite [0-9]+ mass']);
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = expTransform(1, 'xtoa');
end