function [bound, elapsed_time] = demSpmgpJuraBatchKL(file, options, ...
    numActive, iters, numFolds, modelFull)

% DEMSPMGPJURABATCH Demonstrate sparse convolution models on JURA data.

% MULTIGP

dataSetName = ['juraData' file];
display = 0;
bound =  zeros(numFolds,length(numActive));
elapsed_time =  zeros(numFolds,length(numActive));

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
scaleVal = zeros(1,size(yTemp, 2));
biasVal = zeros(1,size(yTemp, 2));
for k =1:size(yTemp, 2),
    biasVal(k) = mean(yTemp{k});
    scaleVal(k) = sqrt(var(yTemp{k}));
end

options.bias = biasVal;
options.scale = scaleVal;

X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));

for i = 1:size(yTemp, 2)
    y{i} = yTemp{i};
    X{i} = XTemp{i};
end

q = 2;
d = size(yTemp, 2);

for seudoPoints =1: length(numActive)
    if seudoPoints == length(numActive),
        sizes = cellfun('size', X, 1);
        options.numActive = sum(sizes);
    else
        options.numActive = numActive(seudoPoints);
    end
    if strcmp(options.approx, 'dtc')
        if options.flagVIKs
            options.nVIKs = options.numActive;
        end
    end
    for fold =1: numFolds
        rand('twister',5^(fold));
        fprintf('Creating MODEL %d DATASET %s APPROX %s NUMACTIVE %d \n',...
        fold, file, options.approx, options.numActive);
        model = multigpCreate(q, d, X, y, options);
        % Initialize parameters
        % Change initial conditions
        modelFullAux = modelFull;
        modelFullAux.paramGroups = speye(size(modelFull.paramGroups,1));
        modelAux = model;
        modelAux.paramGroups = speye(size(model.paramGroups,1));
        [paramsFull, namesFull] = modelExtractParam(modelFullAux);
        indexNoTransform =  paramNameRegularExpressionLookup(modelAux, ' .* sensitivity');
        indexToAugment = paramNameRegularExpressionLookup(modelFullAux, ['multi ' ...
            num2str(length(model.kern.comp)) ' white .* variance']);
        count = length(model.fix);
        count2 = 1;
        for j=1:length(namesFull)           
            if any(j == indexToAugment)
                count2 = count2 + 1;
                nameWhite = ['multi ' num2str(length(model.kern.comp)) ' white ' ... 
                    num2str(count2) ' variance'];
                paramInd = paramNameRegularExpressionLookup(modelAux, nameWhite);            
            else
                paramInd = paramNameRegularExpressionLookup(modelAux, namesFull{j});
            end
            if ~isempty(paramInd)
                count = count + 1;
                if any(paramInd == indexNoTransform)
                    model.fix(count).index = paramInd;
                    model.fix(count).value = paramsFull(j);
                else
                    model.fix(count).index = paramInd;
                    model.fix(count).value = expTransform(exp(paramsFull(j)), 'xtoa');
                end
            end
        end
        % Optimization procedure
        fprintf('Optimizing MODEL %d DATASET %s APPROX %s NUMACTIVE %d \n',...
            fold, file, options.approx, options.numActive);
        tic
        model = multigpOptimise(model, display, iters);       
        elapsed_time(fold, seudoPoints) = toc;
        fprintf('Elapsed time optimization: %f hours\n', elapsed_time(fold,seudoPoints)/3600)
        % Bound        
        bound(fold, seudoPoints) = modelLogLikelihood(model);
        fprintf('Model likelihood: %f\n', bound(fold, seudoPoints));
    end
end




