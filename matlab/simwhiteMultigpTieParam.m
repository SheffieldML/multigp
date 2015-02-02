function tieInd = simwhiteMultigpTieParam(model, options)

% SIMWHITEMULTIGPTIEPARAM Tie parameters for a multigp with SIMWHITE kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses SIMWHITE kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP

for i = 1:options.nlf
    tieInd{i} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* variance']);
end
for i = 1:model.nout
    tieInd{end+1} = paramNameRegularExpressionLookup(model, ['.* ' num2str(i)  ' decay']);
end