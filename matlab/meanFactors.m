function factors = meanFactors(meanFunction, factorType)

% MEANFACTORS Extract factors associated with transformed optimisation space.

% MULTIGP

factors.index = [];
factors.val = [];
if ~isempty(meanFunction.transforms)
  fhandle = str2func([meanFunction.type 'MeanExtractParam']);
  params = fhandle(meanFunction);
  for i = 1:length(meanFunction.transforms)
    index = meanFunction.transforms(i).index;
    factors.index = [factors.index index];
    fhandle = str2func([meanFunction.transforms(i).type 'Transform']);
    factors.val = [factors.val  ...
        fhandle(params(index), factorType)];
  end
end
