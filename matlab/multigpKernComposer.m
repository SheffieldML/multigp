function kernType = multigpKernComposer(type, numOut, numLatent, approx, ...
    specLatent, options)

% MULTIGPKERNCOMPOSER Composes kernel types from some generic options.
% FORMAT
% DESC composes a kernel type to pass to kernel create from a generic
% type and a few options.
% RETURN kernType : an array of cells which indicates the types of kernel
% ARG type : kernel type
% ARG numOut : number of outputs.
% ARG numLatent : number of latent functions.
% ARG specLatent : indicating number for the latent function.
% ARG options : options for kernel. If introduced it's assumed to be
% 'parametric'.
%
% COPYRIGHT : Mauricio Alvarez, 2008
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009
%
% MODIFICATIONS : Mauricio Alvarez, 2010
%
% SEEALSO : multigpCreate

% MULTIGP

if nargin <6
    options = [];
    if nargin<5,
        specLatent = [];
    end
end

switch approx
    case 'ftc'
        kernType = cell(1, numOut+1);
        limit = numOut-numLatent;
    case {'dtc','fitc', 'pitc', 'dtcvar'}
        kernType = cell(1, numOut+numLatent+1);
        limit = numOut;
end

kernType{1} = 'multi';
for i =1:numLatent,
    switch type
        case 'white'
            if strcmp(approx, 'ftc')
                kernType{i+1} = 'none';
            else
                kernType{i+1} = 'white';
                %                kernType{i+1} = 'none';
            end
        otherwise
            kernType{i+1} = 'none';
    end
end


switch type
    case 'gg' % The diffusion kernel.
        if isfield(options, 'kern') && ~isempty(options.kern)
            kernType{specLatent+1} = {'parametric', options.kern, 'gaussian'};
        else
            kernType{specLatent+1} = 'gaussian';
        end
        
        for i = 1:limit
            if isfield(options, 'kern') && ~isempty(options.kern)
                kernType{i+numLatent+1} = {'parametric', options.kern, 'gg'};
            else
                kernType{i+numLatent+1} = 'gg';
            end
        end

    case 'lfm' % the 2nd order differential equation.

        if isfield(options, 'kern') && ~isempty(options.kern)
            kernType{specLatent+1} = {'parametric', options.kern, 'rbf'};
        else
            kernType{specLatent+1} = 'rbf';
        end

        for i = 1:limit
            if isfield(options, 'kern') && ~isempty(options.kern)
                kernType{i+numLatent+1} = {'parametric', options.kern, 'lfm'};
            else
                kernType{i+numLatent+1} = 'lfm';
            end
        end

    case 'lfmwhite' % the 2nd order differential equation with white noise as input

        if strcmp(approx, 'ftc')
            kernType{specLatent+1} = 'white';
        else
            kernType{specLatent+1} = 'rbfinfwhite';
        end
        for i = 1:limit
            kernType{i+numLatent+1} = 'lfmwhite';
        end


    case 'sim' % The 1st order differential equation.

        if isfield(options, 'kern') && ~isempty(options.kern)
            kernType{specLatent+1} = {'parametric', options.kern, 'rbf'};
        else
            kernType{specLatent+1} = 'rbf';
        end

        for i = 1:limit
            if isfield(options, 'kern') && ~isempty(options.kern)
                kernType{i+numLatent+1} = {'parametric', options.kern, 'sim'};
            else
                kernType{i+numLatent+1} = 'sim';
            end
        end

    case 'simwhite' % the 1st order differential equation with white noise as input

        if strcmp(approx, 'ftc')
            kernType{specLatent+1} = 'white';
        else
            kernType{specLatent+1} = 'rbfinfwhite';
        end
        for i = 1:limit
            kernType{i+numLatent+1} = 'simwhite';
        end

    case 'ggwhite'

        if strcmp(approx, 'ftc')
            kernType{specLatent+1} = 'white';
        else
            kernType{specLatent+1} = 'gaussianwhite';
        end

        for i = 1:limit
            kernType{i+numLatent+1} = 'ggwhite';
        end

    case 'rbf' % An independent set of RBF kernels.

        for i = 1:limit
            kernType{i+numLatent+1} = 'ggwhite';
            %kernType{i+numLatent+1} = 'bias';
        end

    case 'white' % Different noise for each output.

        for i = 1:limit
            kernType{i+numLatent+1} = 'white';
        end

    case 'lmc' % The linear model of coregionalization kernel.

        kernType{specLatent+1}= {'parametric', options.kern, 'gaussian'};
        kernType{numLatent+2} = {'parametric', options.kern, 'lmc'};
        kernType = kernType(1:numLatent+2);
        
    case 'whiteblock' % Noise in the LMC kernel
        
        kernType{numLatent+2} = {'parametric', options.kern, 'whiteblock'};
        kernType = kernType(1:numLatent+2);
        
    case 'heat' % The kernel derived from the heat equation
        
        if isfield(options, 'kern') && ~isempty(options.kern)
            if length(options.kern)>1
                kernType{specLatent+1} = {'parametric', options.kern(end), 'rbfh'};
            else
                kernType{specLatent+1} = {'parametric', options.kern, 'rbfh'};
            end
        else
            kernType{specLatent+1} = 'rbfh';
        end

        for i = 1:limit
            if isfield(options, 'kern') && ~isempty(options.kern)
                if length(options.kern)>1
                    kernType{i+numLatent+1} = {'parametric', options.kern(specLatent), 'heat'};
                else
                    kernType{i+numLatent+1} = {'parametric', options.kern, 'heat'};
                end
            else
                kernType{i+numLatent+1} = 'heat';
            end
        end
        
end



