function model = ggwhiteMultigpKernOptions(model, options)

% GGWHITEMULTIGPKERNOPTIONS Changes default options for GGWHITE kernels 
% FORMAT
% DESC Changes default options for the GGWHITE kernel and GAUSSIANWHITE
% kernels
% RETURN model   : model with kernels modified
% ARG    model   : model created
% ARG    options : options for particular kernel
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP
%if strcmp(model.approx,'ftc') || strcmp(model.approx,'dtcvar')
    if isfield(options, 'isArd') && ~isempty(options.isArd)
        optArd  = options.isArd;
    else
        optArd  = false;
    end
    if isfield(options, 'nVIKs') && ~isempty(options.nVIKs)
        if numel(options.nVIKs) ~= numel(options.numActive)
            if numel(options.nVIKs) == 1,
                options.nVIKs = options.nVIKs*ones(1, options.nlf);
            else
                options.nVIKs = options.nVIKs(1)*ones(1, options.nlf);
                warning(['The number of latent functions does not match'...
                    ' the number of elements provided in options.nVIKs']);
            end
        end
        optVIKs = options.nVIKs;  
    else
        optVIKs = ones(1, options.nlf);  
    end
%else
    %error('This kernel should only be used for FTC or DTCVAR approximations')
%end

fhandle1 = str2func([model.kern.comp{1}.comp{1}.type 'KernParamInit']);
fhandle2 = str2func([model.kern.comp{1}.comp{model.nlf+1}.type 'KernParamInit']);
nParamsKernInit = model.kern.nParams;
nParamsInit = 0;
nParamsKern = 0;
for k=1:model.nlf
    nParams = 0;
    if ~strcmp(model.approx,'ftc')        
        model.kern.comp{k}.comp{k} = fhandle1(...
            model.kern.comp{k}.comp{k}, optArd, optVIKs(k));
    end
    nParams = nParams + model.kern.comp{k}.comp{k}.nParams;
    for j=1:model.nout,
        model.kern.comp{k}.comp{model.nlf+j} = fhandle2(model.kern.comp{k}.comp{model.nlf+j}, optArd);
        nParams = nParams + model.kern.comp{k}.comp{model.nlf+j}.nParams;
    end
    nParamsInit = nParamsInit + model.kern.comp{k}.nParams;
    model.kern.comp{k}.nParams = nParams;
    model.kern.comp{k}.paramGroups = speye(nParams);
    nParamsKern = nParamsKern + nParams;
end
model.kern.nParams = nParamsKernInit - nParamsInit + nParamsKern;
model.kern.paramGroups = speye(nParamsKernInit - nParamsInit + nParamsKern);