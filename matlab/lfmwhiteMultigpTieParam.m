function tieInd = lfmwhiteMultigpTieParam(model, options)

% LFMWHITEMULTIGPTIEPARAM Tie parameters for a multigp with LFMWHITE kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses LFMWHITE kernel.
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
    %             tieInd{end+1} = paramNameRegularExpressionLookup(model, ['multi ' ...
    %                 '[0-9]+ lfm ' num2str(i)  ' spring']);
    tieInd{end+1} = paramNameRegularExpressionLookup(model, ['. lfmwhite ' num2str(i)  ' spring']);
    tieInd{end+1} = paramNameRegularExpressionLookup(model, ['multi ' ...
        '[0-9]+ lfmwhite ' num2str(i)  ' damper']);
end