function model = simwhiteMultigpFixParam(model, options)

% SIMWHITEMULTIGPFIXPARAM Fix parameters for a multigp with SIMWHITE kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses SIMWHITE kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP

% This code fixes the latent variances to (1/2)^i (i=0, 1, 2, ..., nlf-1).
count = 0;
for i = 1:model.nlf
    index = paramNameRegularExpressionLookup(model, ...
        ['multi [' num2str(i) '] .* variance']);
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(2^(-i+1), 'xtoa');
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