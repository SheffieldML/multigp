% BATCHHIGHINPUTSGGWHITE Computes DTC VAR bound and full gp likelihood 

% MULTIGP 

clc
clear

rand('twister',1e6);
randn('state',1e6);

indexOut = 2;

inputDim = [2 5 10];
nOuts = [1 2 4 8];
folds = 1; 
nTrainingPoints = 100;
iters = 2000; 

llFull = zeros(length(inputDim),1);
llAprox = cell(length(inputDim),1);

limitNumActive{1} = [20 50 100];
limitNumActive{2} = [80 50 100 200];
limitNumActive{3} = [20 50 100 200 400];
limitNumActive{4} = [20 50 100 200 400 800];

for i = 1:length(inputDim)
    [llFull(i), X, y] = generateDataset('ggwhite', inputDim(i), nOuts(indexOut), nTrainingPoints);
    for j = 1:length(limitNumActive{indexOut}) 
        for k =1:folds
            rand('twister',10^(k));
            randn('state',10^(k));
            llAprox{i}(j,k) =  trainModelHighInGgwhite(X, y, limitNumActive{indexOut}(j), ...
                nOuts(indexOut), iters);
            save(['demBatchHighInputsGgwhite' num2str(nOuts(indexOut))] ,'llFull','llAprox')
        end
    end
end








