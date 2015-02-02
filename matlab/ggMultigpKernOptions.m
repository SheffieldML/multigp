function model = ggMultigpKernOptions(model, options)

% GGMULTIGPKERNOPTIONS Changes the default options for GG kernels 
% FORMAT
% DESC Changes default options for the GG kernel and GAUSSIAN kernels 
% RETURN model   : model with kernels modified
% ARG    model   : model created
% ARG    options : options for particular kernel
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP

if isfield(options, 'isArd') && ~isempty(options.isArd)
    fhandle1 = str2func([model.kern.comp{1}.comp{1}.type 'KernParamInit']);
    fhandle2 = str2func([model.kern.comp{1}.comp{model.nlf+1}.type 'KernParamInit']);
    nParamsKernInit = model.kern.nParams;
    nParamsInit = 0;
    nParamsKern = 0;
    for k=1:model.nlf
        nParams = 0;
        model.kern.comp{k}.comp{k} = fhandle1(...
            model.kern.comp{k}.comp{k}, options.isArd);
        nParams = nParams + model.kern.comp{k}.comp{k}.nParams;
        for j=1:model.nout,
            model.kern.comp{k}.comp{model.nlf+j} = fhandle2(...
                model.kern.comp{k}.comp{model.nlf+j}, options.isArd);
            nParams = nParams + model.kern.comp{k}.comp{model.nlf+j}.nParams;
        end
        nParamsInit = nParamsInit + model.kern.comp{k}.nParams;
        model.kern.comp{k}.nParams = nParams;
        model.kern.comp{k}.paramGroups = speye(nParams);
        nParamsKern = nParamsKern + nParams;
    end
    model.kern.nParams = nParamsKernInit - nParamsInit + nParamsKern;
    model.kern.paramGroups = speye(nParamsKernInit - nParamsInit + nParamsKern);
end
