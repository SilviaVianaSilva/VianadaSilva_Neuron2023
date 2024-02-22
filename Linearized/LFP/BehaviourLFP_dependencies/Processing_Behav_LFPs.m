% Make sure to run these first:

%fullfile(codepath, '/MousePath', 'defineMaze_manual_MH.m')
% Rotate and scale path/maze + add velocity to indata
%fullfile(codepath, '/MousePath', 'runVelocity.m')
%preprocessPathData_Mouse_MH.m


%% Initialize
clc
clear
addpath(genpath(pwd)); % need to modify this more specific
%addpath('C:\Users\Silvia\Desktop\Silvia project code\Tools');
Silvia_PC_path
addpath(genpath(fullfile(codepath,'AnalysisPerTrial','RunAnalysis')));   % e.g. C:\Users\Silvia\Dropbox\SilviaProjectCode

%%  Step 1: Runs just for process_data
% Extract LFP and match to tracking(indata) for S1&S2
%addpath(fullfile(codepath, '/AnalysisPerTrial/'));

process_data = [1450:1457] 
  % perform speed threshold for noramlization of behav to sleep, and extract half-second-binned data for regression analysis 

run(fullfile(codepath,'AnalysisPerTrial','ExtractRegionalEEG_forSpectgAnalysis_MH'))

 
% uncommented here for running sleep first at the same time sept/6/2019

clear
Silvia_PC_path
rmpath(genpath(fullfile(codepath,'AnalysisPerTrial','RunAnalysis')));


%%  Step 2: runs for all animals 
addpath(fullfile(codepath, 'Tools')); %for shaded error bars


for age_gr = 2 %1:2
    if age_gr == 1
        ageism = 'young';
    elseif age_gr == 2
        ageism = 'old';
    end
    
    for rec_area = 1:2% 1:2
        name_exprlayer = 'CA3';       
        array_layers = {'CA1','CA3','DG','DG_CA3'};
        if rec_area ==1
            recordinglayers = [1];
            name_reclayer = 'CA1';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
        elseif rec_area == 2
            recordinglayers = [2,3,4];
            name_reclayer = 'DGCA3';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
        end
[group_A, group_B, animals_analysed, animal_no] = datadef2groups([1,2], recordinglayers, ageism);  %import which animals to analyse, this ignores the process_data above!!
run(fullfile(codepath,'AnalysisPerTrial','ConvertRegionalEEGs.m'))  %normalizes to S1 newer than other function

    end
end


%% Step 3: Cumulate the regional speed in behaviour and output in csv files
%run(fullfile(codepath,'AnalysisPerTrial','RegionalSpeedAccum.m'))  %Runs through all the data inside the script

