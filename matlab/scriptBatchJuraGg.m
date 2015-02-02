clear
clc

file = {'Cd', 'Co', 'Cu', 'Pb'};
numActive =[10 50 100 200 500];
approx = 'pitc';
nlf = 1;
numFolds = 10;
experimentNo = 2;
iters = 1000;
kernType = 'gg';
initialPosition = 'randomComplete';


[maerror, elapsed_time] = demSpmgpJuraBatch(file, numActive, approx, nlf, ...
    numFolds, experimentNo, iters, kernType, initialPosition);