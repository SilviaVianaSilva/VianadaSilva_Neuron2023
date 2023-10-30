% LFP Analysis Sleep
% Extracting LFP during sleep
%
%---------------------------------------------------------
% Matthias Haberl, 2018
%---------------------------------------------------------


disp('Starting LFP Analysis Sleep')

clear
process_data = [405:410]

timer_prog = clock;


%% Animals to process
for iii = process_data    
    clearvars -except process_data iii timer_prog codepath
    clc
    read_PC_path %read in paths
    load(fullfile(codepath, '+Figure8DataOrganization','sessionInfo.mat'));
    All_sessInfo = sessInfo; clear sessInfo
%% Set wavelet parameters
opt.wavelet.dsrate = 16; %Used for spectogramm, 32kHz recording downsampling to 2kHz after wavelet is processed
opt.wavelet.smoothingWin = 20;opt.wavelet.freqRes = 20;
    
    if ~isempty(All_sessInfo(iii).animal)
        
            fprintf('\n -- Starting to process i number: %d --\n', iii);
            close all
            clearvars -except iii process_data All_sessInfo timer_prog aaa animal_numbers opt codepath            
            %%
            sessInfo = All_sessInfo(iii);
            timer_session = tic;
            blocks = {sessInfo.sleepDirs};
            %display(['About to do session ' sessInfo(i).mainDir]);
            for block = sessInfo.sleepDirs
                clearvars power_spect SPG eeg active_SPG quiet_SPG 
                clearvars idx*
                blockDir = fullfile(sessInfo.mainDir, block{1});
                disp(blockDir);fprintf('--> Analyzing: %s\n', block{1});
                
                if sessInfo.age>10
                    Ageism = 'old';
                elseif sessInfo.age<10
                    Ageism = 'young';
                end
                if sessInfo.genotype >=3
                    expressionlayer = 'CA1';
                elseif sessInfo.genotype < 3
                     expressionlayer = 'CA3';
                end
                analysis_dir = ['C:\Users\Silvia\Documents\MATLAB\Sleep_MatFiles_PowerOverFreq\SleepPowerFreq_', Ageism,expressionlayer];mkdir(analysis_dir);
                powfrq_analysisfile = fullfile(analysis_dir,sprintf('animal_%d_i%d_%s.mat',sessInfo.animal,iii, block{1}))
                delete(powfrq_analysisfile);
                exportDir = fullfile(blockDir,'LFP');mkdir(exportDir);
                fullspectfile = fullfile(exportDir,sprintf('Spectogram_channel%d.mat',sessInfo.cellLayerChann));
                fprintf('Clearing %s \n', fullspectfile);
                delete(fullspectfile);
                
                %% Read velocity
                load(fullfile(sessInfo.mainDir,'processedData','indataS.mat'));  
                       
                %% Start processing EEGs
                channelInLayer = sprintf('CSC%d.ncs', sessInfo.cellLayerChann); % Picks the channel that is in the layer
                [eeg, sFreq,~] = readCRTsd(fullfile(blockDir, channelInLayer));
                eeg_raw = Data(eeg);
                eeg_ts = Range(eeg) * 1e-4;
                clearvars eeg
                %% Match EEG and video tracking
                if strcmp(block, 's1'), sleepsess = 1; end
                if strcmp(block, 's2'), sleepsess = 2; end
                TrackingVelocity = indata(sleepsess).v;
                TrackingVelocity = TrackingVelocity';
                TrackingTime = indata(sleepsess).t; 
                clearvars indata
                [~, LFP_idx_nr_of_Vid_start] = (min(abs(TrackingTime(1) - (eeg_ts))));
                [~, LFP_index_nr_of_Video_end] = (min(abs(TrackingTime(end) - (eeg_ts))));
                
                eeg_raw = eeg_raw(LFP_idx_nr_of_Vid_start:LFP_index_nr_of_Video_end);
                eeg_ts = eeg_ts(LFP_idx_nr_of_Vid_start:LFP_index_nr_of_Video_end);
                TrackingTime = TrackingTime - TrackingTime(1);
                TrackingTime = TrackingTime';                

                %% Run Wavelet
                frequs = [2:300];  
                
                fprintf('Recording duration: %s min\n', num2str((eeg_ts(end)-eeg_ts(1))/60));
                Total_time= eeg_ts(end)-eeg_ts(1);  
                if ((eeg_ts(end)-eeg_ts(1))/60)>21;
                errordlg('Warning EEG signal longer than 20min, breaking in two to prevent out-of-memory error');
                
                break_idx=(numel(eeg_ts))/2; 
                disp('Running Wavelet Part 1...');
                [SPG1, t1, f, bandSpecgramFun] = specgramwwd(eeg_raw(1:break_idx,1),sFreq, 2, 300,opt.wavelet);  
                disp('Running Wavelet Part 2...');
                [SPG2, t2, f, bandSpecgramFun] = specgramwwd(eeg_raw(break_idx+1:end,1),sFreq, 2, 300,opt.wavelet);  
                SPG = cat(2,SPG1,SPG2);
                t = [t1, t1(end)+t2];
                clear SPG1 SPG2
                else
                disp('Running Wavelet...')
                [SPG, t, f, bandSpecgramFun] = specgramwwd(eeg_raw,sFreq, 2, 300,opt.wavelet);                       
                end
                
                clearvars eeg_raw
                %% Upsampling of video tracking to match LFP
                TrackingVelocity = imresize(TrackingVelocity,[1,round(size(SPG,2))],'bilinear');
                TrackingTime = imresize(TrackingTime,[1,round(size(SPG,2))],'bilinear');                
                
                exportDir = fullfile(blockDir,'LFP');mkdir(exportDir);
                fullspectfile = fullfile(exportDir,sprintf('Spectogram_channel%d.mat',sessInfo.cellLayerChann))
                save(fullspectfile,'sessInfo','block','SPG','f','t','TrackingVelocity','TrackingTime');
            end
    fprintf('Processing time for Session No.=%d: ',iii);
    toc(timer_session);         
    end
end