function [totalError, elapsed_time_train, elapsed_time_test] = ...
    compilerSparseCoreDTCVAR(dataSetName, options, numTraining, ...
    numActive, display, iters, totFolds)
  
% COMPILERSPARSECOREDTCVAR description.
% MULTIGP

[XTemp, yTemp] = mapLoadData(dataSetName);

load('all_samples.mat');

X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));
XTest = cell(1,size(yTemp, 2));
yTest = cell(1,size(yTemp, 2));
q = size(XTemp{1}, 2);
d = size(yTemp, 2);

totalError = cell(length(numTraining),1);
elapsed_time_train = cell(length(numTraining),1);
elapsed_time_test = cell(length(numTraining),1);

for ns = 1:length(numTraining)
    nSamplesTraining = numTraining(ns);
    totalErrorSamples = zeros(totFolds, size(yTemp, 2), length(numActive{ns}));
    elapsed_time_samples = zeros(totFolds, length(numActive{ns}));
    elapsed_time_samples_test = zeros(totFolds, length(numActive{ns}));
    for j =1:length(numActive{ns}),
        fprintf('NUMBER OF TRAINING POINTS: %d \n',numTraining(ns));
        fprintf('NUMBER OF ACTIVE POINTS: %d \n',numActive{ns}(j));
        options.numActive = numActive{ns}(j);
        options.nVIKs = options.numActive;
        for nFolds =1:totFolds
            fprintf('FOLD NUMBER: %d \n', nFolds);
            for i = 1:size(yTemp, 2)
                y{i} = yTemp{i}(all_samples(1:nSamplesTraining,nFolds),:);
                X{i} = XTemp{i}(all_samples(1:nSamplesTraining,nFolds),:);
                yTest{i} = yTemp{i};
                yTest{i}(all_samples(1:nSamplesTraining,nFolds),:) = [];
                XTest{i} = XTemp{i};
                XTest{i}(all_samples(1:nSamplesTraining,nFolds),:) = [];
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
            [totalErrorSamples(nFolds,:, j), ...
                elapsed_time_samples_test(nFolds, j)] = compilerResults(model,...
                XTest, yTest);
        end
        fprintf('TOTAL TIME: %f \n',sum(elapsed_time_samples_test(:, j)) ...
            + sum(elapsed_time_samples(:, j)));
    end
    totalError{ns} = totalErrorSamples;
    elapsed_time_train{ns} = elapsed_time_samples;
    elapsed_time_test{ns} = elapsed_time_samples_test;
    save(['compiler' upper(options.kernType(1)) options.kernType(2:end) upper(options.approx) 'SeveralVIK']) 
end
%
