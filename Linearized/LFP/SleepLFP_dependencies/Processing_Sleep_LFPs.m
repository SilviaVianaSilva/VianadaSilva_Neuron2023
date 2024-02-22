
%% Make sure to run these first:
% defineMaze_manual.m
% Rotate and scale path/maze + add velocity to indata
% runVelocity.m
% preprocessPathData_Mouse_MH.m


%% Initialize
clc
clear
addpath(genpath(pwd)); % need to modify this more specific

codepath = (pwd);
datapath=('W:\Silvia\Analysis');

%% Extract LFP and match to tracking(indata) for S1&S2
%addpath(fullfile(codepath, '/AnalysisPerTrial/')); 
process_data = [1450:1457] 

addpath(genpath('C:\Users\VianadaSilva_Neuron2023\Linearized\LFP\SleepLFP_dependencies'));

sleep.LFP_Analysis_Sleep_MH  % wavelet on raw LFP to get spectogram
addpath(genpath('C:\Users\Silvia\Dropbox\SilviaProjectCode'));

%process_data = [];
process_data = [1450:1457] %1600:1622];%[3300:3316 3350:3375 3400:3409 3450:3479 3500:3509 3550:3571 3600:3623 3650:3665] 
sleep.process_sleep_spect  % perform speed threshold for noramlization of behav to sleep, and extract half-second-binned data for regression analysis 
rmpath(genpath('C:\Users\Silvia\Dropbox\SilviaProjectCode'));

%