clear
clc

file = {'Cd', 'Co'};
nlf = 1;
numFolds = 10;
experimentNo = 1;
iters = 1000;
kernType = 'ggwhite';

[maerror, elapsed_time] = demJuraBatch(file, nlf, numFolds, experimentNo,...
    iters, kernType);