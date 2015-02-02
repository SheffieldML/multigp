% DEMFOURWALKSWITHLFM Samples from the posterior distribution of a LFM 
% model trained with four sequences of walking and plays the movie for the 
% generated movement.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

clc
clear
close all

dataSetName = 'fourWalks';
experimentNo = 1;
forces = '2';
approx =  'dtcvar';
kernel = 'lfm';
kernel(1) = upper(kernel(1));

capName = dataSetName;
capName(1) = upper(capName(1));
load(['dem' kernel capName num2str(experimentNo) 'Forces' forces 'Approx' approx  '.mat']);

multimodel = model;

indxModelLocal = 1;
model = multimodel.comp{indxModelLocal};
nSamples = 140;

fsp = 120/4;
XTestTemplate = (1:nSamples)'/fsp;
XTest = cell(model.nout+model.nlf,1);
for i=1:model.nout+model.nlf
    XTest{i} = XTestTemplate;
end

% Posterior for the latent functions
kernTemp = model.kern;
kernTemp.approx = 'dtc';
[void, KX_star_X2, Ku_star_u] = globalKernCompute(kernTemp, XTest, {XTest(1:model.nlf), model.X(1:model.nlf)});
mu = cell(1, model.nlf);
varsig = cell(1, model.nlf);
for i=1:model.nlf
    mu{i} = Ku_star_u{i}*model.AinvKuyDinvy{i};
    varsig{i,i} = Ku_star_u{i}*model.Ainv{i,i}*Ku_star_u{i}';
    for j=1:i-1
        varsig{i,j} = Ku_star_u{i}*model.Ainv{i,j}*Ku_star_u{j}';
        varsig{j,i} = varsig{i,j}';
    end
end

samplePosterior = false;
if samplePosterior
    latMu = cell2mat(mu');
    latCovar = cell2mat(varsig);
    latentPos = mat2cell(real(gsamp(latMu, latCovar, 1))', nSamples*ones(model.nlf,1),1);    
else
    latentPos = mu;    
end
[void, KX_star_X2, Ku_star_u] = globalKernCompute(kernTemp, XTest, XTest(1:model.nlf), model.gamma);

KuuinvU = cell(model.nlf,1);
Kuuinv = cell(model.nlf,1);
for i=1:model.nlf 
    Kuuinv{i} = pdinv(Ku_star_u{i}); 
    KuuinvU{i} = Kuuinv{i}*latentPos{i};
end

yP = cell(1, model.nout);
Y = zeros(nSamples, model.nout);
for i=1:model.nout
    yP{i} = zeros(nSamples,1);
    for j=1:model.nlf
        yP{i} = yP{i} + KX_star_X2{i,j}*KuuinvU{j};        
    end
    yP{i} =  yP{i}*model.scale(i) + model.bias(i);
    Y(:,i) = yP{i}; 
end
    
channels = Y;
load('skelForWalkingWithLfm.mat');
%skel = acclaimReadSkel('../../datasets/matlab/mocap/cmu/12/12.asf');
%[tmpchan, skel] = acclaimLoadChannels('../../datasets/matlab/mocap/cmu/12/12_01.amc', skel);

figure(1)
set(gcf, 'Position', [ 67 277 560 420]);
limits = [-5 6; -5 5; -1 29];
walkSamplePlayData(skel, channels, limits, 1/50)
%walkSamplePlayData(skel, channels, [], 1/200)

