function model = simwhiteMultigpKernOptions(model, options)

% SIMWHITEMULTIGPKERNOPTIONS Changes the default options for SIMWHITE kernels 
% FORMAT
% DESC Changes default options for the SIMWHITE kernel and RBF kernels in the
% multigp structure.
% RETURN model   : model with kernels modified
% ARG    model   : model created
% ARG    options : options for particular kernel
%
% COPYRIGHT : David Luengo, Mauricio A. Alvarez, 2010

% MULTIGP


if isfield(options, 'positiveTime') && ~isempty(options.positiveTime)
    for i=1:model.nlf
        for j=1:model.nlf
            model.kern.comp{i}.comp{j}.positiveTime = options.positiveTime;
        end        
        for j=1:model.nout
            model.kern.comp{i}.comp{model.nlf + j}.positiveTime = options.positiveTime;
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
