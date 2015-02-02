% DEMTOYCOMPARISONICMRANK2 Generates samples from an ICM with Rank 2
% FORMAT
% DESC Generates samples from a multi-output Gaussian process with
% covariance given by the intrinsic coregionalization model. This script
% was used to generate figure 4.2, in the paper "Kernels for vector-valued
% functions: a review", appearing in Foundations and Trends in Machine
% Learning.
%
% COPYRIGHT : Mauricio A. Alvarez, 2012

% MULTIGP


clc
clear
close all
rand('twister', 2196620689); 
randn('seed', 1.424690759000000e+09);

N = 500;
x = linspace(0,5, N)';
nSamples = 2;
D = 2;
rankC = 2;
nlf = 1;
kernType = cell(nlf+1,1);

kernType{1} = 'cmpnd';
for i=1:nlf
   kernType{i+1} = {'parametric', struct('nout', D, 'rankCorregMatrix', rankC), 'lmc'}; 
end

kern = kernCreate(x, kernType);
B = [1 0.5; 0.5 1.5];

lengthScale = 0.5;
precisionU = 1./(lengthScale.^2);

for i=1:nlf
   kern.comp{i}.B = B;
   kern.comp{i}.precisionU = precisionU(i); 
end

KICM = kernCompute(kern, x);
yICM = gsamp(zeros(D*N,1), KICM, nSamples);
yI_ICM = cell(D,nSamples);

for ns = 1:nSamples
    startVal = 1;
    endVal = 0;
    for d =1:D
        endVal = endVal + N;
        yI_ICM{d,ns} = (yICM(ns, startVal:endVal))';        
        startVal = endVal + 1;
    end
end

lWidth = 2;
spec = {'-k', '--k'};
fontsize = 15;

figure
for ns = 1:nSamples
    subplot(2,1,1)
    hold on
    a1 = plot(x, yI_ICM{1,ns}, spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('ICM $R_q=2$, $f_1(x)$', 'interpreter', 'latex', 'fontsize', 1.2*fontsize)
    subplot(2,1,2)
    hold on
    a3 = plot(x, yI_ICM{2,ns},spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('ICM $R_q=2$, $f_2(x)$', 'interpreter', 'latex',  'fontsize', 1.2*fontsize)
end

%print('-depsc', '~/Dropbox/MOLSurvey/FoTML/toyICMPRank2', '-loose');





