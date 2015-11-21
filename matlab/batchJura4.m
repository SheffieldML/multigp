% BATCHJURA4  Batch of an Independent Gaussian Process using the Jura Dataset. 

% MULTIGP 

% In this demo, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.
% When changing the kernel, the fix values in this code and the indexes in
% the tieParam vector in multigpCreate must also be changed.


clc
clear
% n1 = rand('seed');
% n2 = randn('seed');
% %
file = {'Cd', 'Co', 'Cu', 'Pb'};
dataSetName = 'data_jura';
experimentNo = 1;
ntrainX =200;
ntrainX2 =50;
iters = 1000;
approx = 'none';
saveFigures=1;
% load data



% Inclusion of the latent function
options.isSparse = 0;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 1;  % Indicates if the outputs are to be standarized
options.outKernelName = 'gg'; % Indicates the name of the output kernel: lfm, sim or gg

total_t = zeros(10,4); 
maerror = zeros(10,4);

for q =1:4
    for k = 1:10,
        fprintf('Experiment: %s Cross-Validation: %f\n',file{q},k);
        rand('seed',(q+k+2)*10^6);
        randn('seed',(q+k+2)*10^6);
        data = mapLoadData([dataSetName '_' file{q}]);
        data.nin = 1;  % Number of latent functions
        data.nout = 1;
        data.Xtrain = data.Xtrain(1);
        data.Ytrain = data.Ytrain(:,1);
        data.Xtest = data.Xtest(1);  
        data.Ytest = data.Ytest(:,1);
        missingData = cell(data.nout,1);
        [total_t(k,q), maerror(k,q)] = batch_demJura(options, data, ntrainX, ntrainX2, approx, missingData, iters);
        save('./results/Jura4Ind','total_t','maerror');
    end
end



