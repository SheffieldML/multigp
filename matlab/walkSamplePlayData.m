function walkSamplePlayData(skelStruct, channels, limits, frameLength) 

% WALKSAMPLEPLAYDATA Play sampled walking motion.
% FORMAT 
% DESC plays channels from a motion capture skeleton and channels.
% ARG skelStruct : the skeleton for the motion.
% ARG channels : the channels for the motion.
% ARG limits : limits to plot the axes
% ARG frameLength : the framelength for the motion.
%
% COPYRIGHT : Mauricio A. Alvarez, Neil D. Lawrence, 2010
%
% SEEALSO : skelPlayData, acclaimPlayData

% MULTIGP

if nargin < 4
  frameLength = 1/120;
end

clf
handle = skelVisualise(channels(1, :), skelStruct);

xlim = [limits(1,1) limits(1,2)];
ylim = [limits(2,1) limits(2,2)];
zlim = [limits(3,1) limits(3,2)];
set(gca, 'xlim', xlim, ...
         'ylim', ylim, ...
         'zlim', zlim);
title('Walking motion', 'FontSize', 15);

% Play the motion
for j = 1:size(channels, 1)
  pause(frameLength)
  skelModify(handle, channels(j, :), skelStruct);
end
