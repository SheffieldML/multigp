function tieInd = lfmMultigpTieParam(model, options)

% LFMMULTIGPTIEPARAM Tie the parameters for a multigp model with LFM kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses LFM kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP

for i = 1:options.nlf
    tieInd{i} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* inverse width']);
end
for i = 1:model.nout
    tieInd{end+1} = paramNameRegularExpressionLookup(model, ['. lfm ' num2str(i)  ' spring']);
    tieInd{end+1} = paramNameRegularExpressionLookup(model, ['multi ' ...
        '[0-9]+ lfm ' num2str(i)  ' damper']);
end