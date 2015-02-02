function options = multigpOptions(approx)

% MULTIGPOPTIONS Return default options for the MOCAP examples in the LFM model.
% FORMAT
% DESC returns the default options in a structure for a MULTIGP model.
% ARG approx : approximation type, either 'none' (no approximation),
% 'fitc' (fully
% independent training conditional) or 'pitc' (partially
% independent training conditional.
% RETURN options : structure containing the default options for the
% given approximation type.
%
% SEEALSO : multigpCreate
%
% COPYRIGHT : Mauricio Alvarez, 2008
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MULTIGP

  options.type = 'multigp';
  
  if nargin<1
    options.approx = 'ftc';
  else
    options.approx = approx;
  end
  % Basic kernel
  options.kernType = 'gg';
  % Include noise in the model
  options.includeNoise = true;
  %Learn the scales
  options.learnScales = false;
  % Include independent kernel in the model
  options.includeInd = false;
  % Include options to tie the parameters
  options.tieOptions = multigpTieOptions;
  % Method for optimization
  options.optimiser = 'scg';
  % One latent function.
  options.nlf = 1;
  
  % Set to a given mean function to have a mean function.
  options.meanFunction = [];
  % Options structure for mean function options.
  options.meanFunctionOptions = [];
  
  
  
  switch options.approx
   case 'ftc'
    options.numActive = [];
   case {'dtc','fitc','pitc', 'dtcvar'}
    options.numActive = 15;
    options.fixInducing = true;
    options.fixIndices = 1:options.numActive;
    options.includeScalesp = 0;
    options.tieInducing = false;
    if (strcmp(options.approx, 'dtc') || strcmp(options.approx, 'dtcvar'))
        options.beta = 1e-3;
    else
        options.beta = 1e3;
    end
    % Initial position of the inducing variables. Options are 'random',
    % in which random initial locations taken from the used data;
    % 'espaced' the initial locations are equally spaced chosen across
    % all dimensions; 'fixIndices' the indices in the fixIndices
    % options are employed; 'kmeans' the kmeans method gives the
    % initial positions.        
    options.initialInducingPositionMethod = 'random'; 
    if strcmp(options.approx, 'dtcvar')
        % Learns the sensitivities variationally
        options.varS = false;
    end        
  end
end


function options = multigpTieOptions
  options.tieIndices = false;
  options.selectMethod = 'free';
end

