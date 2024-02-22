disp('Starting LFP Analysis Sleep')

failed_iii = [];
timer_prog = clock;
%for aaa = animal_numbers
% find

%% Animals to process
for iii = process_data    
    clearvars -except process_data iii failed_iii timer_prog codepath
    clc
    Silvia_PC_path %read in paths
    load(fullfile(codepath, '+Figure8DataOrganization','sessionInfo.mat'));
    All_sessInfo = sessInfo; clear sessInfo
%% Set wavelet paraameters
opt.wavelet.dsrate = 16; %Used for spectogramm, 32kHz recording downsampling to 2kHz after wavelet is processed
opt.wavelet.smoothingWin = 20;opt.wavelet.freqRes = 20;
    
    if ~isempty(All_sessInfo(iii).animal)
        
        %try
            fprintf('\n -- Starting to process i number: %d --\n', iii);
            close all
            clearvars -except iii process_data All_sessInfo timer_prog aaa animal_numbers failed_iii opt codepath            
            %%
            sessInfo = All_sessInfo(iii);
            timer_session = tic;
            blocks = {sessInfo.sleepDirs};
            %display(['About to do session ' sessInfo(i).mainDir]);
            for block = sessInfo.sleepDirs
                clearvars power_spect SPG eeg active_SPG quiet_SPG 
                clearvars idx*
                % try
                blockDir = fullfile(sessInfo.mainDir, block{1});
                disp(blockDir);fprintf('--> Analyzing: %s\n', block{1});
                
                %% define infos, for filenames; and delete old files first! otherwise errors can be masked later on, when re-running analysis 
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
                eeg_raw = Data(eeg);%eeg_raw = eeg_raw(1:1000000,1);
                eeg_ts = Range(eeg) * 1e-4;
                clearvars eeg
                %% Fix mismatch of EEG and video tracking
                if strcmp(block, 's1'), sleepsess = 1; end
                if strcmp(block, 's2'), sleepsess = 2; end
                TrackingVelocity = indata(sleepsess).v;
                TrackingVelocity = TrackingVelocity';
                TrackingTime = indata(sleepsess).t; %- indata(sleepsess).t(1); %set to zero   
                clearvars indata
                %[eeg_ts, eeg_raw] = eeg.AlignLFPwInData(TrackingTime,eeg_ts, eeg_raw)
                [~, LFP_idx_nr_of_Vid_start] = (min(abs(TrackingTime(1) - (eeg_ts))));
                [~, LFP_index_nr_of_Video_end] = (min(abs(TrackingTime(end) - (eeg_ts))));
                
                if (LFP_idx_nr_of_Vid_start==LFP_index_nr_of_Video_end), errordlg(sprintf('Mismatch of Video tracking and EEG data timestamps in s%s in i:%s, video-ts = %i and eeg-ts = %i ?!?', num2str(sleepsess), num2str(iii), TrackingTime(1), eeg_ts(1)));
                 warndlg(sprintf('Failed days: %s', num2str(failed_iii)), 'Summary failed days');
                    continue
                end
                eeg_raw = eeg_raw(LFP_idx_nr_of_Vid_start:LFP_index_nr_of_Video_end);
                eeg_ts = eeg_ts(LFP_idx_nr_of_Vid_start:LFP_index_nr_of_Video_end);
                TrackingTime = TrackingTime - TrackingTime(1);
                TrackingTime = TrackingTime';                
                        %{
                        %% Extract first minute EEG to save
                        [~, idx1]=min(abs(eeg_ts-(eeg_ts(1)+60)));
                        eeg_firstMin  = eeg_raw(1:idx1,1);
                        eeg_ts_firstMin  = eeg_ts(1:idx1,1);
                        %% Extract last 1min of sleep for this
                        [~, idx]=min(abs(eeg_ts-(eeg_ts(end)-60)));
                        eeg_lastMin  = eeg_raw(idx:end,1);
                        eeg_ts_lastMin  = eeg_ts(idx:end,1);
                        %}                                             
                        %{
                        %% Downsampling
                        eeg_ds = downsample(eeg_raw,4);
                        eeg_ts_ds = downsample(eeg_ts,4);
                        sFreq_ds = sFreq/4;
                        %}                
                %% Run Wavelet
                frequs = [2:300];  
                
                fprintf('Recording duration: %s min\n', num2str((eeg_ts(end)-eeg_ts(1))/60));
                Total_time= eeg_ts(end)-eeg_ts(1);  
                if ((eeg_ts(end)-eeg_ts(1))/60)>21;
                errordlg('Warning EEG signal longer than 20min, breaking in two');
                
                break_idx=(numel(eeg_ts))/2; 
                disp('Running Wavelet Part 1...');
                [SPG1, t1, f, bandSpecgramFun] = specgramwwd(eeg_raw(1:break_idx,1),sFreq, 2, 300,opt.wavelet);  %changed now to use all eeg_raw instead of just last minute
                disp('Running Wavelet Part 2...');
                [SPG2, t2, f, bandSpecgramFun] = specgramwwd(eeg_raw(break_idx+1:end,1),sFreq, 2, 300,opt.wavelet);  %changed now to use all eeg_raw instead of just last minute               
                SPG = cat(2,SPG1,SPG2);
                t = [t1, t1(end)+t2];
                clear SPG1 SPG2
                else
                disp('Running Wavelet...')
                [SPG, t, f, bandSpecgramFun] = specgramwwd(eeg_raw,sFreq, 2, 300,opt.wavelet);  %changed now to use all eeg_raw instead of just last minute                                 
                end
                
                clearvars eeg_raw
                %% Upsampling of video tracking to match LFP
                TrackingVelocity = imresize(TrackingVelocity,[1,round(size(SPG,2))],'bilinear');
                TrackingTime = imresize(TrackingTime,[1,round(size(SPG,2))],'bilinear');                
                     
                
                exportDir = fullfile(blockDir,'LFP');mkdir(exportDir);
                fullspectfile = fullfile(exportDir,sprintf('Spectogram_channel%d.mat',sessInfo.cellLayerChann))
                save(fullspectfile,'sessInfo','block','SPG','f','t','TrackingVelocity','TrackingTime');%,'quiet_power_spect', 'active_power_spect', 'quiet_power_std', 'active_power_std','active_t','quiet_t','total_t')
            end
    fprintf('Processing time for Session No.=%d: ',iii);
    toc(timer_session);         
    end
    
if ~isempty(failed_iii);     warndlg(sprintf('Failed days: %s', num2str(failed_iii)), 'Summary failed days');    end    
end