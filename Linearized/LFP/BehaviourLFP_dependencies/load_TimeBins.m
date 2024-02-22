% load_TimeBins
%Loads opt.binninByTime
% Bin size is used for run_LFP_Analysis_inParsedMaze_MH
% and later for plots

load(fullfile('C:\Users\Silvia\Dropbox\SilviaProjectCode\+Figure8DataOrganization','TimeSpentInPasedMaze.mat'),'TimeSpentInMaze','AverageTime');

opt.binninByTime.d0 = round(50 * [AverageTime.d0])% Total time = 15.1282);  %inserted here average time in seconds per area on Fig8, divided by sum of the time
%opt.binninByTime.d0 = 1000/opt.binninByTime.TotalnumberBins * [286,  82, 161, 104, 367];
opt.binninByTime.d2 = round(50 * [AverageTime.d2]);% Total time= 18.3352);  %inserted here average time in seconds per area on Fig8, divided by sum of the time
opt.binninByTime.d10 =  round(50 * [AverageTime.d10]);% Total time = /25.5904);

