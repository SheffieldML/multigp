% BATCHJURA3  Batch of the Sparse Multi Output Gaussian Process using the PITC approx over the Jura datatset. 

% MULTIGP

%In this demo, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.
% When changing the kernel, the fix values in this code and the indexes in
% the tieParam vector in multigpCreate must also be changed.

clc
clear
n1 = rand('seed');
n2 = randn('seed');
%
% Add necessary toolboxes
file = {'Cd', 'Co', 'Cu', 'Pb'};
dataSetName = 'data_jura';
experimentNo = 1;
ntrainX =200;
ntrainX2 =[10 50 100 200 500];
iters = 30;
approx = 'pitc';
saveFigures=1;
% load data



% Inclusion of the latent function
options.isSparse = 1;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 1;  % Indicates if the outputs are to be standarized
options.outKernelName = 'gg'; % Indicates the name of the output kernel: lfm, sim or gg
total_t = zeros(10,4,length(ntrainX2)); 
maerror = zeros(10,4,length(ntrainX2));

for r=1:length(ntrainX2),
    for q =1:4,
        for k = 1:10,
            data = mapLoadData([dataSetName '_' file{q}]);
            data.nin = 1;  % Number of latent functions
            missingData = cell(data.nout,1);
            [total_t(k,q,r), maerror(k,q,r)] = batch_demJura(options, data, ntrainX, ntrainX2(r), approx, missingData, iters);
            save('./results/Jura3','total_t','maerror', 'n1', 'n2');
        end
    end
end
