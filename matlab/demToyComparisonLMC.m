% DEMTOYCOMPARISONLMC Samples from a multi GP with LMC covariance
% FORMAT
% DESC Generates samples from a multi-output Gaussian process with a
% covariance function that follows the Linear Model of Coregionalization 
% (LMC), with rank 2 and 2 latent functions. This script was used to 
% generate figure 4.4, in the paper "Kernels for vector-valued functions: 
% a review", appearing in Foundations and Trends in Machine Learning.
%
% COPYRIGHT : Mauricio A. Alvarez, 2012

% MULTIGP

clc
clear
close all
rand('twister', 3406384884); 
randn('seed', 578868733);

N = 500;
x = linspace(0,5, N)';
nSamples = 2;
D = 2;
rankC = 2;
nlf = 2;
kernType = cell(nlf+1,1);

kernType{1} = 'cmpnd';
for i=1:nlf
   kernType{i+1} = {'parametric', struct('nout', D, 'rankCorregMatrix', rankC), 'lmc'}; 
end

kern = kernCreate(x, kernType);

lengthScales = [1 0.2];
precisionU = 1./(lengthScales.^2);

A{1} = [0.1 0.5;0.2 0.5];
A{2} = [1 0.3;0.5 0.2]; 

B{1} = [1.4 0.5; 0.5 1.2];

B{2} = [1 0.5; 0.5 1.3];


for i=1:nlf
   kern.comp{i}.A = A{i};
   kern.comp{i}.B = B{i};
   kern.comp{i}.precisionU = precisionU(i); 
end

KLMC = kernCompute(kern, x);
yLMC = gsamp(zeros(D*N,1), KLMC, nSamples);
yI_LMC = cell(D,nSamples);


for ns = 1:nSamples
    startVal = 1;
    endVal = 0;
    for d =1:D
        endVal = endVal + N;
        yI_LMC{d,ns} = (yLMC(ns, startVal:endVal))';        
        startVal = endVal + 1;
    end
end

lWidth = 2;
spec = {'-k', '--k'};
fontsize = 15;
figure

means = [-1.5 2.5; -1.5 2.2];

for ns = 1:nSamples
    subplot(2,1,1)
    hold on
    a1 = plot(x, yI_LMC{1,ns} - mean(yI_LMC{1,ns}) + means(1, ns), spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('LMC with $R_q=2$ and $Q=2$, $f_1(x)$', 'interpreter', 'latex', 'fontsize', 1.2*fontsize)
    subplot(2,1,2)
    hold on
    a3 = plot(x, yI_LMC{2,ns} - mean(yI_LMC{2,ns}) + means(2, ns),spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('LMC with $R_q=2$ and $Q=2$, $f_2(x)$', 'interpreter', 'latex',  'fontsize', 1.2*fontsize)
end

%print('-depsc', '~/Dropbox/MOLSurvey/FoTML/toyLMC', '-loose');





