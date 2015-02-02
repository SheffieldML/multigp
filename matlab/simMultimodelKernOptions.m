function model = simMultimodelKernOptions(model, options)

% SIMMULTIMODELKERNOPTIONS Changes the default options for SIM kernels
% FORMAT
% DESC Changes default options for the SIM kernel and RBF kernels in the
% spmultimodel structure.
% RETURN model   : model with kernels modified
% ARG    model   : model created
% ARG    options : options for particular kernel
%
% COPYRIGHT : David Luengo, Mauricio A. Alvarez, 2009

% MULTIGP


if ~model.varS
    if isfield(options, 'isNormalised') && ~isempty(options.isNormalised)
        for i=1:model.nlf
            model.kern.comp{i}.comp{1}.isNormalised = options.isNormalised;
            for j=1:model.numOutGivenLat(i)
                model.kern.comp{i}.comp{1 + j}.isNormalised = options.isNormalised;
            end
        end
    end
    if isfield(options, 'isStationary') && ~isempty(options.isStationary)
        for i=1:model.nlf
            model.kern.comp{i}.comp{1}.isStationary = options.isStationary;
            for j=1:model.numOutGivenLat(i)
                model.kern.comp{i}.comp{1 + j}.isStationary = options.isStationary;
            end
        end
    end
    if isfield(options, 'isNegativeS') && ~isempty(options.isNegativeS)
        for i=1:model.nlf
            for j=1:model.numOutGivenLat(i)
                model.kern.comp{i}.comp{1 + j}.isNegativeS = options.isNegativeS;
                if options.isNegativeS == true
                    model.kern.comp{i}.comp{1 + j}.transforms.index = [1 2];
                else
                    model.kern.comp{i}.comp{1 + j}.transforms.index = [1 2 3];
                end
            end
        end
    end
else
    if isfield(options, 'isNormalised') && ~isempty(options.isNormalised)
        for i=1:model.nlf
            model.kern.comp{i}.comp{1}.isNormalised = options.isNormalised;
            for j=1:model.nout
                model.kern.comp{i}.comp{1 + j}.isNormalised = options.isNormalised;
            end
        end
    end
    if isfield(options, 'isStationary') && ~isempty(options.isStationary)
        for i=1:model.nlf
            model.kern.comp{i}.comp{1}.isStationary = options.isStationary;
            for j=1:model.nout
                model.kern.comp{i}.comp{1 + j}.isStationary = options.isStationary;
            end
        end
    end
    if isfield(options, 'isNegativeS') && ~isempty(options.isNegativeS)
        for i=1:model.nlf
            for j=1:model.nout
                model.kern.comp{i}.comp{1 + j}.isNegativeS = options.isNegativeS;
                if options.isNegativeS == true
                    model.kern.comp{i}.comp{1 + j}.transforms.index = [1 2];
                else
                    model.kern.comp{i}.comp{1 + j}.transforms.index = [1 2 3];
                end
            end
        end
    end
end

