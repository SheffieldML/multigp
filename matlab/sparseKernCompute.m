function model = sparseKernCompute(model)

% SPARSEKERNCOMPUTE Computes the kernels for the sparse approximation in the convolution process framework.
% FORMAT
% DESC computes the kernels for the sparse approximation in the convolution process framework.
% This  version is more computational effcient, because for the sparse
% approximations we don't need to evaluate any cross covariances in the ouputs.
% In a set up with many outputs, this should be a waste of resources if we use
% the calssical way of computing in the kerntoolbox.
% On the other hand, we can decide how to store the covariances of the
% outputs: as a vector, corresponding to the diagonal of the
% covariances of the outputs for the fitc or as a full covariance,
% corresponding to the pitc.
% ARG model : the model structure
% RETURN model : the modified model structure with the kernels updated.
%
% COPYRIGHT : Mauricio Alvarez, 2008, 2010

% MULTIGP

% Compute the Kuu and Kyu parts

if strcmp(model.kernType, 'lmc')    
    for r = 1:model.nlf,
        model.Kuu{r,1} = multiKernComputeBlock(model.kern.comp{r},  model.X{r}, r, r);
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            model.KuuGamma{r,1} = model.Kuu{r,1} + model.gamma(r)*eye(size(model.Kuu{r,1}));
        end
        for i=1:model.nout
            model.Kyu{i,r} = real(multiKernComputeBlock(model.kern.comp{r}, ...
                model.X{i+model.nlf}, model.X{r}, r, r));
            switch model.approx
                case {'fitc', 'dtcvar'}
                    model.Kyy{i,r} = kernDiagCompute(model.kern.comp{r}.comp{r}, model.X{i+model.nlf});
                case 'pitc'
                    model.Kyy{i,r} = multiKernComputeBlock(model.kern.comp{r}, model.X{i+model.nlf}, r, r);
            end
        end
    end
else
    for r = 1:model.nlf,
        model.Kuu{r,1} = zeros(model.k(r));
        for c = 1: length(model.kern.comp)
            model.Kuu{r,1} = model.Kuu{r,1} + real(multiKernComputeBlock(model.kern.comp{c},  model.X{r}, r, r));
        end
        for i =1:model.nout,
            model.Kyu{i,r} = real(multiKernComputeBlock(model.kern.comp{r},  model.X{i+model.nlf},...
                model.X{r}, i+model.nlf, r));
        end
    end
    % Compute the Kyy part
    for j =1:model.nout,
        switch model.approx
            case {'fitc', 'dtcvar'}
                model.Kyy{j,1} = zeros(size(model.X{model.nlf + j},1),1);
            case 'pitc'
                model.Kyy{j,1} = zeros(size(model.X{model.nlf + j},1));
        end
        for c = 1: length(model.kern.comp)
            switch model.approx
                case {'fitc', 'dtcvar'}
                    model.Kyy{j,1} = model.Kyy{j,1} + real(kernDiagCompute(model.kern.comp{c}.comp{j+model.nlf}, ...
                        model.X{j+model.nlf}));
                case 'pitc'
                    model.Kyy{j,1} = model.Kyy{j,1} + real(multiKernComputeBlock(model.kern.comp{c}, ...
                        model.X{j+model.nlf}, j+model.nlf, j+model.nlf));
            end
        end
    end
end
