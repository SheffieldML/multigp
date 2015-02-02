% DEMTOYCOMPARISONSLFM Samples from a multi GP with LMC covariance
% FORMAT
% DESC Generates samples from a multi-output Gaussian process with a
% covariance function that follows the Linear Model of Coregionalization 
% (LMC), with rank 1 and 2 latent functions. This covariance corresponds to
% the Semiparametric Latent Factor Model covariance (SLFM). The script was 
% used to generate figure 4.3, in the paper "Kernels for vector-valued 
% functions: a review", appearing in Foundations and Trends in Machine 
% Learning.
%
% COPYRIGHT : Mauricio A. Alvarez, 2012

% MULTIGP


clc
clear
close all
rand('twister', 120676318); 
randn('seed', 1.648668321000000e+09);

N = 500;
x = linspace(0,5, N)';
nSamples = 2;
D = 2;
rankC = 1;
nlf = 2;
kernType = cell(nlf+1,1);

kernType{1} = 'cmpnd';
for i=1:nlf
   kernType{i+1} = {'parametric', struct('nout', D, 'rankCorregMatrix', rankC), 'lmc'}; 
end

kern = kernCreate(x, kernType);

lengthScales = [1 0.2];

precisionU = 1./(lengthScales.^2);
A{1} = [0.5; 1];
A{2} = [1;0.5]; 

for i=1:nlf
   kern.comp{i}.A = A{i};
   kern.comp{i}.B = kern.comp{i}.A*kern.comp{i}.A';
   kern.comp{i}.precisionU = precisionU(i); 
end


KSLFM = kernCompute(kern, x);
ySLFM = gsamp(zeros(D*N,1), KSLFM, nSamples);
yI_SLFM = cell(D,nSamples);


for ns = 1:nSamples
    startVal = 1;
    endVal = 0;
    for d =1:D
        endVal = endVal + N;
        yI_SLFM{d,ns} = (ySLFM(ns, startVal:endVal))';        
        startVal = endVal + 1;
    end
end

lWidth = 2;
spec = {'-k', '--k'};
fontsize = 15;
figure
ameans = [0 2];
for ns = 1:nSamples
    subplot(2,1,1)
    hold on
    a1 = plot(x, yI_SLFM{1,ns} - mean(yI_SLFM{1,ns}) + ameans(ns), spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('LMC with $R_q=1$ and $Q=2$, $f_1(x)$', 'interpreter', 'latex', 'fontsize', 1.2*fontsize)
    subplot(2,1,2)
    hold on
    a3 = plot(x, yI_SLFM{2,ns} - mean(yI_SLFM{2,ns})  + ameans(ns),spec{ns}, 'linewidth', lWidth);
    set(gca, 'fontname', 'times', 'fontsize', fontsize , 'color', 'none');
    xlabel('LMC with $R_q=1$ and $Q=2$, $f_2(x)$', 'interpreter', 'latex',  'fontsize', 1.2*fontsize)
end

%print('-depsc', '~/Dropbox/MOLSurvey/FoTML/toySLFM', '-loose');





