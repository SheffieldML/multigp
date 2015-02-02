% DEMAISTATS Reproduces demo presented at AISTATS 2009

% MULTIGP

colordef white

disp('Ready ... play training data.')
r = input('Type ''R'' to run or ''S'' to skip: ', 's');
switch r
  case {'r', 'R'}
   close all
   clear all
   nReps = 1;
   rep = 1;
   nSamples = 100;

   while(rep)
     load 'data18'     
     ind = round(linspace(1, size(channels, 1), nSamples));
     figure(1)
     set(gcf, 'Position', [ 67 277 560 420]);
     
     for j=1:nReps
       clf
       title('Motion 18 (train)')
       balancePlayData(skel, channels(ind, :), limits, 18, 49, 1/25);
       pause(1.5)
     end
     figure(2)
     set(gcf, 'Position', [671 272 560 420]);
     load 'data19'
     ind = round(linspace(1, size(channels, 1), nSamples));
          
     for j=1:nReps
       clf
       title('Motion 19 (train)')
       balancePlayData(skel, channels(ind, :), limits, 19, 49, 1/25);
       pause(1.5)
     end
     r2 = input('Type ''R'' to repeat or ''C'' to continue: ', 's');
     switch r2
      case {'r', 'R'}
       rep = 1;
      otherwise
       rep = 0;
     end
   end
 otherwise
end

disp('Ready ... play test data.')
r = input('Type ''R'' to run or ''S'' to skip: ', 's');
switch r
  case {'r', 'R'}
   close all
   clear all
   nReps = 1;
   rep = 1;
   nSamples = 100;

   while(rep)
     load 'data18'
     ind = round(linspace(1, size(channels, 1), nSamples));
     figure(3)
     set(gcf, 'Position', [671 272 560 420]);
     

     for j=1:nReps
       clf
       title('Motion 20 (test)')
       balancePlayData(skel, channels(ind, :), limits, 20, 49, 1/25);
       pause(1.5)
     end
     r2 = input('Type ''R'' to repeat or ''C'' to continue: ', 's');
     switch r2
      case {'r', 'R'}
       rep = 1;
      otherwise
       rep = 0;
     end
   end
 otherwise
end





disp('Ready ... simulate data.')
r = input('Type ''R'' to run or ''S'' to skip: ', 's');
switch r
  case {'r', 'R'}
   close all
   clear all

   load additionalCmu49BalanceArmMod2.mat
   
   lfmResultsDynamic('LfmCmu49BalanceArm', 27, 'balance', initPos, skel, muMod, muMod2);
 otherwise
end
