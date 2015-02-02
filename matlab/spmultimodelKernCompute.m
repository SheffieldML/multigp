function model = spmultimodelKernCompute(model)

% SPMULTIMODELKERNCOMPUTE Computes the kernels in the sparse multimodel
% FORMAT
% DESC Computes the kernels for the sparse multimodel.
% ARG model : model with previous kernels. 
% RETURN model : model with updated kernels.
%
% COPYRIGHT : Mauricio Alvarez, 2010

% MULTIGP


if model.varS
    switch model.isSpeed
        case 1
            basett=  cputime;
            if model.subSpeed == 1
                fhandle = str2func([model.kernType 'KernCompute']);
                if isfield(model, 'gamma') && ~isempty(model.gamma)
                    [model.Kyy, model.Kyu, model.KuuGamma] = fhandle(model.kern, ...
                        model.outX, model.latX, model.gamma);
                else
                    [model.Kyy, model.Kyu, model.Kuu] = fhandle(model.kern, ...
                        model.outX, model.latX);
                end
            else
                fhandle = str2func([model.kernType 'KernCompute2']);
                [model.Kyy, model.Kyu, model.KuuGamma] = fhandle(model.kern, ...
                    model.outX, model.latX, model.sizeX, model.k, model.gamma);
            end
            endtt = cputime - basett;
            a = 1;
         case 2
            % Compute Kuu
            for k = 1:model.nlf
                model.Kuu{k,1} = real(model.kernFuncNames.computeLat(model.kern.comp{k}.comp{1},  model.latX{k}));
                if isfield(model, 'gamma') && ~isempty(model.gamma)
                    model.KuuGamma{k,1} = model.Kuu{k,1} + model.gamma(k)*eye(size(model.Kuu{k,1}));
                end
            end
            for i = 1:model.nout,
                for j = 1: model.nlf
                    % Compute Kff
                    model.Kyy{i,j} = real(model.kernFuncNames.computeOut(model.kern.comp{j}.comp{1+i}, ...
                        model.outX{i}));
                    % Compute Kfu, which corresponds to K_{\hat{f}}u, really.
                    model.Kyu{i,j} = real(model.kernFuncNames.computeCross(model.kern.comp{j}.comp{1+i}, ...
                        model.kern.comp{j}.comp{1}, model.outX{i},...
                        model.latX{j}));
                end
            end
        case 3
            % Compute Kuu
            for k = 1:model.nlf
                model.Kuu{k,1} = real(multiKernComputeBlock(model.kern.comp{k},  model.latX{k}, 1, 1));
                if isfield(model, 'gamma') && ~isempty(model.gamma)
                    model.KuuGamma{k,1} = model.Kuu{k,1} + model.gamma(k)*eye(size(model.Kuu{k,1}));
                end
            end
            for i = 1:model.nout,
                for j = 1: model.nlf
                    % Compute Kff
                    model.Kyy{i,j} = real(kernDiagCompute(model.kern.comp{j}.comp{1+i}, ...
                        model.outX{i}));
                    % Compute Kfu, which corresponds to K_{\hat{f}}u, really.
                    %model.Kyu{i,j} = model.qs.mu(i,j)*real(multiKernComputeBlock(model.kern.comp{j},  model.outX{i},...
                    %    model.latX{j}, i+1, 1));
                    model.Kyu{i,j} = real(multiKernComputeBlock(model.kern.comp{j},  model.outX{i},...
                        model.latX{j}, i+1, 1));
                end
            end
        otherwise
            error('Unknown SPEED option')
    end
else
    latGivenOut = zeros(1, model.nout);
    for k = 1:model.nlf
        model.Kuu{k,1} = real(multiKernComputeBlock(model.kern.comp{k},  model.latX{k}, 1, 1));
        if isfield(model, 'gamma') && ~isempty(model.gamma)
            model.KuuGamma{k,1} = model.Kuu{k,1} + model.gamma(k)*eye(size(model.Kuu{k,1}));
        end
        %%%
        %model.Kuu{k,1} = model.Kuu{k,1} + 0.01*eye(model.k);
        %%%
        whichOutputs = model.indOutGivenLat{k};
        latGivenOut(whichOutputs) = latGivenOut(whichOutputs) + 1;
        for j =1:model.numOutGivenLat(k);
            temp = real(multiKernComputeBlock(model.kern.comp{k},  model.outX{whichOutputs(j)},...
                model.latX{k}, j+1, 1));
            model.Kuy{k,1}{j,1} = temp';
            model.Kyu{whichOutputs(j),1}{latGivenOut(whichOutputs(j)),1} = temp;
        end
    end
    for i = 1:model.nout,
        latentFunc = model.indLatGivenOut{i};
        switch model.approx
            case {'dtc', 'fitc', 'dtcvar'}
                model.Kyy{i,1} = zeros(model.sizeX(i),1);
            case 'pitc'
                model.Kyy{i,1} = zeros(model.sizeX(i));
        end
        for j = 1: length(latentFunc)
            whichOutputs = model.indOutGivenLat{latentFunc(j)};
            localOutput = find(i == whichOutputs );
            switch model.approx
                case {'dtc','fitc', 'dtcvar'}
                    model.Kyy{i,1} = model.Kyy{i,1} + real(kernDiagCompute(model.kern.comp{latentFunc(j)}.comp{1+localOutput}, ...
                        model.outX{i}));
                case 'pitc'
                    model.Kyy{i,1} = model.Kyy{i,1} + real(multiKernComputeBlock(model.kern.comp{latentFunc(j)}, ...
                        model.outX{i}, 1+localOutput, 1+localOutput));
            end
        end
    end
end