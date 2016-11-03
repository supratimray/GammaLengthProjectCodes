% Readme.m

% How to get MP software
% Matching pursuit software with stochastic dictionary and the supporting
% matlab files are available at http://www.fuw.edu.pl/~durka/software/mp/.
% The file 'mp31.exe' (which runs the matching pursuit algorithm) should be
% added to the working directory. Some minor corrections has to be done in
% original ‘gabor.m’ function provided by Prof.Piotr Durka and colleagues.
% The modified version of ‘gabor.m’ is included in this folder.

% How to generate the plots

% 1. Simply run "plotFigure1" to generate Figure 1.

% 2. Plots in the second figure are generated using the program
% ‘testPerformanceSynthData’, which has three inputs. 
% a. thresholdFractionList: array of threshold factors
% b. electrodeNum: %86-low, 83-medium, 29-high;
% c. numMeanBursts: number of bursts to be injected per trial. If you assign ‘numMeanBursts’ as an empty matrix, the number of bursts is calculated from a Poisson distribution). 
% Note: MP computation is very slow, and can be disabled by setting the variable ‘showMPResults’ on line 32 to 0.      
% For example, testPerformanceSynthData([0.5 0.1 0.25 0.75],83,1) will
% generate Figures 2E, J, K and L.

% 3. run runGetBurstLengthRealData.m to generate Figure 3. Analysis is done
% only on the three electrodes for which data is provided. 