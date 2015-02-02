function model = simMultigpKernOptions(model, options)

% SIMMULTIGPKERNOPTIONS Changes the default options for SIM kernels 
% FORMAT
% DESC Changes default options for the SIM kernel and RBF kernels 
% RETURN model   : model with kernels modified
% ARG    model   : model created
% ARG    options : options for particular kernel
%
% COPYRIGHT : David Luengo, Mauricio A. Alvarez, 2009

% MULTIGP

if isfield(options, 'isNormalised') && ~isempty(options.isNormalised)
    for i=1:model.nlf
        for j=1:model.nlf
            model.kern.comp{i}.comp{j}.isNormalised = options.isNormalised;
        end        
        for j=1:model.nout
            model.kern.comp{i}.comp{model.nlf + j}.isNormalised = options.isNormalised;
        end
    end
end

if isfield(options, 'isStationary') && ~isempty(options.isStationary)    
    for i=1:model.nlf
        for j=1:model.nlf
            model.kern.comp{i}.comp{j}.isStationary = options.isStationary;
        end
        for j=1:model.nout
            model.kern.comp{i}.comp{model.nlf+j}.isStationary = options.isStationary;
        end
    end
end
    
