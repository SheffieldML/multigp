function [yEst, z, W] = demPpcaFxData(dataSetName, nlf);

% DEMPPCAFXDATA Demonstrate PPCA model (for comparison purposes with LFM
% model) on exchange rates data.
%
% FORMAT
% DESC Provides the results of PPCA applied to the FX data set with missing
% data. By default the dimension of the latent space is equal to the number
% of outputs or the rank of W, whichever is lower.
% RETURN yEst : estimated outputs.
% RETURN z : transformed outputs in latent space.
% RETURN W : transformation matrix such that yEst = z * W'.
% ARG dataSetName : The data set to load.
%
% FORMAT
% DESC Does the same as above but allows the specification of the dimension
% of the latent space.
% RETURN yEst : estimated outputs.
% RETURN z : transformed outputs in latent space.
% RETURN W : transformation matrix such that yEst = z * W'.
% ARG dataSetName : The data set to load.
% ARG nlf : Optional parameter. Dimension of the latent space. If not used
% the dimension of the latent space is equal to the number of outputs. If
% nlf is higher than the rank of W, then nlf is set to this value.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, loadFxData, fxDataResults

% MULTIGP


rand('seed', 1e6);
randn('seed', 1e6);

% Load data

data = loadFxData(dataSetName);

% Options for the demo

optionMeanVal = 'sampleMean';
% optionMeanVal = 'initVal';
offsetX = 1;
display = 1;
niter = 1000;

nout = size(data.y, 2);
ndata = size(data.y{1}, 1);

% Set the missing data for the demo

missingData = cell(1, nout);
missingData{4} = 50:100;
% missingData{5} = 200:250;
missingData{6} = 100:150;
% missingData{7} = 1:50;
missingData{9} = 150:200;

% Set the inputs and outputs in the correct format

X = zeros(ndata, nout);
y = zeros(ndata, nout);
meanVal = zeros(1, nout);
scaleVal = ones(1, nout);
count = offsetX + (0:ndata-1)';
for i = 1:nout
    if ~isempty(missingData{i})
        data.y{i}(missingData{i}) = 0.0;
    end
    ind1 = find(data.y{i} ~= 0.0);
    ind2 = find(data.y{i} == 0.0);
    data.y{i}(ind2) = NaN;
    scaleVal(i) = std(data.y{i}(ind1));
    if strcmp(optionMeanVal, 'sampleMean')
        meanVal(i) = mean(data.y{i}(ind1));
    elseif strcmp(optionMeanVal, 'initVal')
        meanVal(i) = data.y{i}(ind1(1));
    else
        error('Mean value option not defined');
    end
    y(:, i) = (data.y{i} - meanVal(i))/scaleVal(i);
    X(:, i) = count;
end

% Embedding of the missing data with PPCA

if (nargin == 2)
    nEmbed = nlf;
else
    nEmbed = nout;
end

[z, sigma2, W] = ppcaEmbed(y, nEmbed);
if rank(W)<nEmbed
    nEmbed = rank(W);
    warning(['Rank of W lower than latent dimension. Latent dimension reduced to ' num2str(nEmbed)]);
    [z, sigma2, W] = ppcaEmbed(y, nEmbed);
end

% Creating and training an independent GP model for each output

options = gpOptions('ftc');

for i=1:size(z, 2)
    disp(' ');
    disp(['%%% Output ' num2str(i) ' %%%']);
    eval(['model' num2str(i) ' = gpCreate(1, 1, X(:, i), z(:, i), options);']);
    eval(['model' num2str(i) ' = gpOptimise(model' num2str(i) ', display, niter);']);
end

% Recover the estimate of the original data

yEst = z*W';
yEst = yEst .* repmat(scaleVal, size(y, 1), 1) + repmat(meanVal, size(y, 1), 1);

if nargin == 2
    eval(['save demPpca' dataSetName 'Nlf' num2str(nlf) ' yEst z W;']);
else
    eval(['save demPpca' dataSetName 'Nlf' num2str(nout) ' yEst z W;']);
end
    