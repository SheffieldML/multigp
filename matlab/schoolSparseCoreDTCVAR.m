function [totalError, elapsed_time_samples] = ...
    schoolSparseCoreDTCVAR(dataSetName, options, numActive, display, iters, totFolds)

% SCHOOLSPARSECOREDTCVAR description.

% MULTIGP

[XTemp, yTemp] = mapLoadData(dataSetName);



X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));
XTest = cell(1,size(yTemp, 2));
yTest = cell(1,size(yTemp, 2));
q = size(XTemp{1}, 2)-1;
d = size(yTemp, 2);

elapsed_time_samples = zeros(totFolds, length(numActive));
totalError = zeros(totFolds, length(numActive));

for j =1:length(numActive),
    fprintf('NUMBER OF ACTIVE POINTS: %d \n',numActive(j));
    options.numActive = numActive(j);
    options.nVIKs = options.numActive;
    rand('twister',10^7);
    randn('state',10^7);            
    for nFolds =1:totFolds
        fprintf('FOLD NUMBER: %d \n', nFolds);        
        for i = 1:size(yTemp, 2)
            maxSamples = size(yTemp{i}, 1);
            randIndex = randperm(maxSamples);
            samplesTraining = floor(0.75*maxSamples);
            indexTraining = randIndex(1:samplesTraining);
            indexTesting = randIndex(samplesTraining+1:end);
            y{i} = yTemp{i}(indexTraining);
            X{i} = XTemp{i}(indexTraining,1:q);
            options.nRepeats{i} = XTemp{i}(indexTraining,q+1);
            yTest{i} = yTemp{i}(indexTesting,:);
            XTest{i} = XTemp{i}(indexTesting,1:q);
            options.bias(i) = mean(y{i});
            options.scale(i) = std(y{i});
        end
        % Creates the model
        model = multigpCreate(q, d, X, y, options);
        %Trains the model
        tic;
        model = multigpOptimise(model, display, iters);
        elapsed_time_samples(nFolds, j) = toc;
        % Save the results.
        CC = schoolResults(model, XTest, yTest);
        totalError(nFolds, j) = mean(CC);
    end
    fprintf('TOTAL TIME: %f \n',sum(elapsed_time_samples(:, j)));
    save(['school' upper(options.kernType(1)) options.kernType(2:end) upper(options.approx) 'SeveralVIK']) 
end
