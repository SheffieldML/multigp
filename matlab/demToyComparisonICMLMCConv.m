% DEMTOYCOMPARISONICMLMCCONV Samples from multi-output covariances
% FORMAT
% DESC Generates samples from a multi-output Gaussian process with
% different covariance functions, namely the intrinsic coregionalization
% model (ICM) with rank 1, the linear model of coregionalization (LMC) with
% rank 2 and 2 latent functions, and the process convolution covariance
% with 1 latent function. This script was used to generate figure 5.1, in 
% the paper "Kernels for vector-valued functions: a review", appearing in 
% Foundations and Trends in Machine Learning.
%
% COPYRIGHT : Mauricio A. Alvarez, 2012

% MULTIGP


clc
clear
close all
rand('twister', 3406384884); 
randn('seed', 1.116547193000000e+09);

N = 500;
x = linspace(0,5, N)';
nSamples = 2;
D = 2;
rankC = 1;

kernGG = kernCreate(x, {'multi', 'ggwhite', 'ggwhite'});
kernICM = kernCreate(x, {'parametric', struct('nout', D, 'rankCorregMatrix', rankC), 'lmc'});
kernICM.A = [1;5];
kernICM.B = kernICM.A*kernICM.A';
kernICM.precisionU = 1/(0.5^2); 

lengthScales = [1 0.1];
precisionG = 1./(lengthScales.^2);
varianceGG = [1 1];

for d =1:D
    kernGG.comp{d}.precisionG = precisionG(d);
    kernGG.comp{d}.variance = varianceGG(d);
end


rankLMC = 2;
nlfLMC = 2;
kernType = cell(nlfLMC+1,1);

kernType{1} = 'cmpnd';
for i=1:nlfLMC
   kernType{i+1} = {'parametric', struct('nout', D, 'rankCorregMatrix', rankLMC), 'lmc'}; 
end

lengthScales = [1 0.2];
kernLMC = kernCreate(x, kernType);
precisionU = 1./(lengthScales.^2);

A{1} = [0.1 0.5;0.2 0.5];
A{2} = [1 0.3;0.5 0.2]; 

B{1} = [1.4 0.5; 0.5 1.2];
B{2} = [1 0.5; 0.5 1.3];

for i=1:nlfLMC
   kernLMC.comp{i}.A = A{i};
   kernLMC.comp{i}.B = B{i};
   kernLMC.comp{i}.precisionU = precisionU(i); 
end

KLMC = kernCompute(kernLMC, x);
KGG = kernCompute(kernGG, x);
KICM = lmcKernCompute(kernICM, x);

yGG = gsamp(zeros(D*N,1), KGG, nSamples);
yICM = gsamp(zeros(D*N,1), KICM, nSamples);

rand('twister', 3406384884); 
randn('seed', 578868733);
yLMC = gsamp(zeros(D*N,1), KLMC, nSamples);

yI_LMC = cell(D,nSamples);
yI_GG = cell(D,nSamples);
yI_ICM = cell(D,nSamples);


for ns = 1:nSamples
    startVal = 1;
    endVal = 0;
    for d =1:D
        endVal = endVal + N;
        yI_GG{d,ns} = (yGG(ns, startVal:endVal))';
        yI_ICM{d,ns} = (yICM(ns, startVal:endVal))';
        yI_LMC{d,ns} = (yLMC(ns, startVal:endVal))';
        startVal = endVal + 1;
    end
end

lWidth = 2;
spec = {'-k', '--k'};
fontsize = 15;

meansGG = [0 1.5; 0 3];
means = [-1.5 2.5; -1.5 2.2];
meansICM = [0 2; 0 5];
xlim = [0 5];

for ns = 1:nSamples
    subplot(2,3,1)
    hold on
    a1 = plot(x, yI_ICM{1,ns} - mean(yI_ICM{1,ns}) + meansICM(1, ns), spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none',  'xlim', xlim);
    xlabel('ICM, $f_1(x)$', 'interpreter', 'latex', 'fontsize', fontsize)
    subplot(2,3,2)
    hold on
    a2 = plot(x, yI_LMC{1,ns} - mean(yI_LMC{1,ns}) + means(1, ns), spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none', 'xlim', xlim);
    xlabel('LMC, $f_1(x)$', 'interpreter', 'latex', 'fontsize', fontsize)
    subplot(2,3,3)
    hold on
    a3 = plot(x, yI_GG{1,ns} - mean(yI_GG{1,ns}) + meansGG(1,ns), spec{ns},'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none', 'xlim', xlim);
    xlabel('PC, $f_1(x)$', 'interpreter', 'latex', 'fontsize', fontsize)
    %%%%%%%%%%%
    subplot(2,3,4)
    hold on
    a4 = plot(x, yI_ICM{2,ns} - mean(yI_ICM{2,ns}) + meansICM(2, ns) ,spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none', 'xlim', xlim);
    xlabel('ICM, $f_2(x)$', 'interpreter', 'latex',  'fontsize', fontsize)    
    subplot(2,3,5)
    hold on
    a5 = plot(x, yI_LMC{2,ns} - mean(yI_LMC{2,ns}) + means(2, ns),spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none', 'xlim', xlim);
    xlabel('LMC, $f_2(x)$', 'interpreter', 'latex',  'fontsize', fontsize)
    subplot(2,3,6)
    hold on
    a6 = plot(x, yI_GG{2,ns} - mean(yI_GG{2,ns}) + meansGG(2,ns), spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none', 'xlim', xlim);
    xlabel('PC, $f_2(x)$', 'interpreter', 'latex',  'fontsize', fontsize)
end

%print('-depsc', '~/Dropbox/MOLSurvey/FoTML/toyICMLMCPCv1', '-loose');





