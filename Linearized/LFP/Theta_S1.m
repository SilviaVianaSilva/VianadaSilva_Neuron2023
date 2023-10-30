function Theta_S1

% Testing for  S1 theta frequency
%---------------------------------------------------------
% Matthias Haberl, 2023
%---------------------------------------------------------


addpath('Metadata');
addpath('Tools');
addpath('ProjectCode\');
addpath(genpath('RunAnalysis')); 
load('.\+Figure8DataOrganization\sessionInfo.mat');
All_sessInfo = sessInfo; clear sessInfo

expressionlayer = 'CA3';

for age_gr = 1:2
    if age_gr == 1
        ageism = 'young';
    elseif age_gr == 2
        ageism = 'old';
    end
    
    for rec_area = 1:2
        name_exprlayer = 'CA3';
        array_layers = {'CA1','CA3','DG','DG_CA3'};
        if rec_area ==1
            recordinglayers = [1];
            name_reclayer = 'CA1';  
        elseif rec_area == 2
            recordinglayers = [2,3,4];
            name_reclayer = 'DGCA3';  
        end       
        [group_A, group_B, animals_analysed, animal_no] = datadef2groups([1,2], recordinglayers, ageism);
        animals_groupA =  unique(animal_no(group_A));
        animals_groupB =  unique(animal_no(group_B));
        animals = {animals_groupA; animals_groupB};
        
        for groupAB = 1:2
            
            all_theta_active_S1 = [];
            cum_animal_numbers = [];
            for aaa = animals{groupAB}  % Animals to process
                fprintf('Analyzing animal number: %s \n', num2str(aaa));
                
                i_list = find((animal_no==aaa) & (group_A | group_B))  % groups generated in datadef2groups, like this recording layers aren't mixed etc.
                is_age = All_sessInfo(i_list(1)).age; %Check the age of this animal
                reclayer = All_sessInfo(i_list(1)).tListLoc;
                
                for iii = i_list

                 %   try
                    fprintf('Analyzing i number: %s \n', num2str(iii));
                    %%  Get LFP
                    fprintf('Getting LFP \n');
                    lfp = loadEEG(All_sessInfo(iii).mainDir, 's1', All_sessInfo(iii).cellLayerChann);
                    
                       %eegfilename = fullfile((All_sessInfo(iii).mainDir), 's1', strcat('CSC',num2str(All_sessInfo(iii).cellLayerChann),'.ncs'));
                      %[tst, samp]=   getRawCSCData(eegfilename,1,1000)
                    %% Get tracking and velocity
                    indata_file = fullfile(All_sessInfo(iii).mainDir, 'processedData','indataS.mat');
                    fprintf('Getting indata %s \n', indata_file);
                    load(indata_file);
                    
                    % TrackingTime = indata.t;  eeg_ts = lfp.ts; eeg_raw =  lfp.samp;
                    %% Align tracking and LFP
                    eeg_ts =  lfp.ts;
                    eeg_raw =  lfp.samp;
                    indata = indata(1);
                    [eeg_ts, eeg_raw] = COX_Project.cox_functions.AlignLFPwInData(indata.t, eeg_ts, eeg_raw);
                    sFreq =  lfp.sampFreq;
                    clearvars lfp
                    
                    %% Run Wavelet
                    opt.wavelet.dsrate = 1;  %Used for spectogramm
                    opt.wavelet.smoothingWin = 20;opt.wavelet.freqRes = 30;
                    frequs = [2:40];
                    
                    fprintf('Recording duration: %s min\n', num2str((eeg_ts(end)-eeg_ts(1))/60));
                    Total_time= eeg_ts(end)-eeg_ts(1);
                    
                    if ((eeg_ts(end)-eeg_ts(1))/60)>21;
                        errordlg('Warning EEG signal longer than 20min, breaking in two');
                        break_idx=(numel(eeg_ts))/2;
                        disp('Running Wavelet Part 1...');
                        [SPG1, t1, f, bandSpecgramFun] = specgramwwd(eeg_raw(1:break_idx,1),sFreq, 2, 40,opt.wavelet);  %changed now to use all eeg_raw instead of just last minute
                        disp('Running Wavelet Part 2...');
                        [SPG2, t2, f, bandSpecgramFun] = specgramwwd(eeg_raw(break_idx+1:end,1),sFreq, 2, 40,opt.wavelet);  %changed now to use all eeg_raw instead of just last minute
                        SPG = cat(2,SPG1,SPG2);
                        t = [t1, t1(end)+t2];
                        clear SPG1 SPG2
                    else
                        disp('Running Wavelet...')
                        [SPG, t, f, bandSpecgramFun] = specgramwwd(eeg_raw,sFreq, 2, 40, opt.wavelet);  %changed now to use all eeg_raw instead of just last minute
                    end

                    clearvars eeg_raw
                    
                   
                    %% Upsampling of video tracking to match LFP
                    TrackingTime = indata.t';
                    TrackingVelocity = indata.v';
                    TrackingVelocity = imresize(TrackingVelocity,[1,round(size(SPG,2))]);
                    TrackingTime = imresize(TrackingTime,[1,round(size(SPG,2))]);
                    
                    
                    

                    
                    %% Brainstate                    
                    disp('Analyzing Brainstate...')
                    Theta_band = [6,12];Delta_band = [1.9,4];
                    theta_freq = ((f>=Theta_band(1))&(f<Theta_band(2)));
                    Theta = (SPG(theta_freq,:));
                    delta_freq = ((f>Delta_band(1))&(f<Delta_band(2)));
                    Delta = (SPG(delta_freq,:));
                    brainstate = (max(Theta,[],1)./max(Delta,[],1));
 
                    %% Scaling to seconds
                    downsampledTrackingVelocity = imresize(TrackingVelocity,[1,round(t(end)/2)]);
                    downsampledTrackingTime = imresize(TrackingTime',[1,round(t(end))]);
                    active_t = sum(downsampledTrackingVelocity>2);
                    quiet_t = sum(downsampledTrackingVelocity<=2);
                    figure
                    L1 = plot(downsampledTrackingVelocity,'b'); hold on
                    L2 = plot(1000*imresize(max(Theta,[],1),[1,round(t(end)/2)]),'r'); hold on
                    L3 = plot(1000*imresize(max(Delta,[],1),[1,round(t(end)/2)]),'g'); hold on
                    L4 = plot(imresize((brainstate),[1,round(t(end)/2)]),'k'); hold off
                    legend([L1,L2,L3,L4],{'Veloc','Theta','Delta','Th/Del'});
                    title(sprintf('Active time: %s sec --- Quiet time: %s sec', num2str(active_t),num2str(quiet_t)));
                    print(fullfile(blockDir,'Brainstate_Velocity_inSec.pdf'), '-dpdf');
                    close
                     
                    %% Seperate SPG by velocity
                    quiet_SPG = SPG(:,(TrackingVelocity<2) ); 
                    active_SPG = [];
                    active_SPG = SPG(:,(TrackingVelocity>=2));
                    quiet_t = (t(3)-t(2)) * size(quiet_SPG,2); %calculate actual amount of time in quiet
                    active_t = (t(3)-t(2)) * size(active_SPG,2); %numel(t(:,(mean(Theta,1)./mean(Delta,1))>2));
                    total_t = (t(3)-t(2)) * size(SPG,2);
                    fprintf('--> Active time: %s sec --- Quiet time: %s sec\n <--', num2str(active_t),num2str(quiet_t));
                    
                    %% Check with
                    power_spect =  mean(SPG,2);
                    quiet_power_spect = mean(quiet_SPG,2);
                    active_power_spect = mean(active_SPG,2);
                    quiet_power_std = std(quiet_SPG,0,2);
                    active_power_std = std(active_SPG,0,2);
                    shadedErrorBar(f,active_power_spect,active_power_std,{'k-','markerfacecolor','k'},1); hold on
                    shadedErrorBar(f,quiet_power_spect,quiet_power_std,{'r-','markerfacecolor','r'},1); hold off

                    set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 16, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                    ylabel('Mean Power'); xlabel('Freq (Hz)');
                    xlim(gca,[2 300]);
                    ylim(gca,[0 0.0025]);
                    print(fullfile(blockDir,'powerfrequ.pdf'), '-dpdf');
                    close
                    
                    [PKS,LOCS] = findpeaks(active_power_spect, 'SortStr','descend');
                    theta_peaks = find(f(LOCS)>6 & f(LOCS)<12);
                    height = PKS(theta_peaks);
                    
                    theta_freq = f(LOCS(theta_peaks(1)));
                    
                    animal_active_theta_freq = cat(2, animal_active_theta_freq, theta_freq);
                    animal_quiet_power_spect = cat(2,animal_quiet_power_spect, quiet_power_spect);
                    animal_active_power_spect = cat(2,animal_active_power_spect,active_power_spect);
                    clear quiet_power_spect active_power_spect quiet_power_std active_power_std
                         
                end
                close all hidden


                frequ_analysis.plot_Behav_quietSl_MovSl(filename,f,mean(animal_spect,3),mean(animal_sleepNormalized,3),mean(animal_quiet_power_spect,2),mean(animal_active_power_spect,3));
                
                outfolder = ('U:\Silvia\Analysis\LFP_Analysis\S1_frequ');mkdir(outfolder);
                filename = fullfile(outfolder,sprintf('Animal_%d',aaa));
                save(filename,'animal_quiet_power_spect','f','animal_active_power_spect');
                
                all_theta_active_S1 = cat(1,all_theta_active_S1, mean(animal_active_theta_freq));
                cum_animal_numbers = cat(1,cum_animal_numbers, aaa);
            end
            
            write_csvfile(outfolder, sprintf('Genotype_%d_%s_%s-recordings_active_S1_theta', groupAB, ageism,name_reclayer ) ,all_theta_active_S1)
            write_csvfile(outfolder, sprintf('Genotype_%d_%s_%s-recordings_animal', groupAB, ageism,name_reclayer ) , cum_animal_numbers)
            
        end
        
        
    end
end

end



function lfp = loadEEG(infolder, task, channel)
eegfilename = fullfile(infolder, task, strcat('CSC',num2str(channel),'.ncs'));
fprintf('Reading: %s\n', eegfilename);
[eeg, sFreq] = readCRTsd(eegfilename);% does the data get flipped?

lfp.samp = Data(eeg)*-1; %flip data;

lfp.samp=downsample(lfp.samp,16);  %downsample; take every 16th data point
lfp.ts = Range(eeg);
lfp.ts = downsample(lfp.ts,16);
lfp.sampFreq = sFreq/16;
end




