% BATCHJURA1 Batch of the Full Multi Output Gaussian Process using the Jura Dataset. 

% MULTIGP

%In this demo, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.
% When changing the kernel, the fix values in this code and the indeces in
% the tieParam vector in multigpCreate must also be changed.


clc
clear
% n1 = rand('seed');
% n2 = randn('seed');
%
% Add necessary toolboxes
% file = {'Cd', 'Co', 'Cu', 'Pb'};
file = {'Cu', 'Co', 'Cd', 'Pb'};
dataSetName = 'data_jura';
experimentNo = 1;
ntrainX =200;
ntrainX2 =50;
iters = 1000;
approx = 'none';
saveFigures=1;




% Inclusion of the latent function
options.isSparse = 0;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 1;  % Indicates if the outputs are to be standarized
options.outKernelName = 'rbf'; % Indicates the name of the output kernel: lfm, sim or gg

total_t = zeros(10,4); 
maerror = zeros(10,4);

for q =1:4
    for k = 1:10,
        data = mapLoadData([dataSetName '_' file{q}]);
        data.nin = 1;  % Number of latent functions
        missingData = cell(data.nout,1);
        [total_t(k,q), maerror(k,q)] = batch_demJura(options, data, ntrainX, ntrainX2, approx, missingData, iters);
        save('./results/Jura4Ind.mat','total_t','maerror');
    end
end



