% DEMTOYCOMPARISONICMRANK1 Generates samples from an ICM with Rank 1
% FORMAT
% DESC Generates samples from a multi-output Gaussian process with
% covariance given by the intrinsic coregionalization model with rank 1. 
%
% COPYRIGHT : Mauricio A. Alvarez, 2012

% MULTIGP

clc
clear
close all
rand('twister', 782593036); 
randn('seed', 1.048404045000000e+09);

N = 500;
x = linspace(0,5, N)';
nSamples = 2;
D = 2;
rankC = 1;
nlf = 1;
kernType = cell(nlf+1,1);

kernType{1} = 'cmpnd';
for i=1:nlf
   kernType{i+1} = {'parametric', struct('nout', D, 'rankCorregMatrix', rankC), 'lmc'}; 
end

kern = kernCreate(x, kernType);

A = [1;5];
lengthScale = 0.5;
precisionU = 1./(lengthScale.^2);

for i=1:nlf
   kern.comp{i}.A = A;
   kern.comp{i}.B = kern.comp{i}.A*kern.comp{i}.A';
   kern.comp{i}.precisionU = precisionU(i); 
end

bkern = kernCreate(x, kern.comp{1}.basicKernelType);
bkern.precisionU = kern.comp{1}.precisionU;
Klatent = kernCompute(bkern, x);

yb = gsamp(zeros(N,1), Klatent, nSamples);
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
    xlabel('ICM $R_q=1$, $f_1(x)$', 'interpreter', 'latex', 'fontsize', 1.2*fontsize)
    subplot(2,1,2)
    hold on
    a3 = plot(x, yI_ICM{2,ns},spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('ICM $R_q=1$, $f_2(x)$', 'interpreter', 'latex',  'fontsize', 1.2*fontsize)
end

%print('-depsc', '~/Dropbox/MOLSurvey/FoTML/toyICMRank1', '-loose');





