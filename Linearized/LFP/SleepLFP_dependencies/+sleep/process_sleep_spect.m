% process_sleep_spect = Step 2
% 03/31/2018 MH

clc
disp('Starting')

failed_iii = [];
timer_prog = clock;

%% Animals to process
for iii = process_data
    clearvars -except process_data iii failed_iii timer_prog
    clc
    Silvia_PC_path %read in paths
    load(fullfile(codepath, '+Figure8DataOrganization','sessionInfo.mat'));
    All_sessInfo = sessInfo; clear sessInfo
    
    if ~isempty(All_sessInfo(iii).animal)
        
        %try
        fprintf('\n -- Starting to process i number: %d --\n', iii);
        close all
        clearvars -except iii process_data All_sessInfo timer_prog aaa animal_numbers failed_iii opt
       
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
            
            %% Load full spectogram
            spectDir = fullfile(blockDir,'LFP');mkdir(spectDir);
            fullspectfile = fullfile(spectDir,sprintf('Spectogram_channel%d.mat',sessInfo.cellLayerChann))
            load(fullspectfile,'sessInfo','block','SPG','f','t','TrackingVelocity','TrackingTime');
            
            
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

            %% For normalizing behaviour to Sleep
            
            %% Seperate SPG by velocity (but not any more by brainstate)
            quiet_SPG = SPG(:,(TrackingVelocity<2) ); % (brainstate<1)& % 31/03/2018: removed brainstate as parameter to account
            active_SPG = [];
            active_SPG = SPG(:,(TrackingVelocity>=2));% (brainstate>2)&
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
            %active=plot(f,mean(active_SPG,2),'g');hold on; quiet=plot(f,mean(quiet_SPG,2),'r'); hold off;
            %shadedErrorBar
            shadedErrorBar(f,active_power_spect,active_power_std,{'k-','markerfacecolor','k'},1); hold on
            shadedErrorBar(f,quiet_power_spect,quiet_power_std,{'r-','markerfacecolor','r'},1); hold off
            %legend({'stddev';'Active';'stddev';'Quiet'});
            set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 16, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
            ylabel('Mean Power'); xlabel('Freq (Hz)');
            xlim(gca,[2 300]);
            ylim(gca,[0 0.0025]);
            print(fullfile(blockDir,'powerfrequ.pdf'), '-dpdf');
            close
           
            analysis_dir = ['W:\Silvia\Analysis\LFP_Analysis\Sleep_MatFiles_PowerOverFreq\SleepPowerFreq_', Ageism,expressionlayer];mkdir(analysis_dir);
            powfrq_analysisfile = fullfile(analysis_dir,sprintf('animal_%d_i%d_%s.mat',sessInfo.animal,iii, block{1}))
            delete(powfrq_analysisfile);
            save(powfrq_analysisfile, 'power_spect','quiet_power_spect','active_power_spect','quiet_power_std', 'active_power_std','f','sessInfo','block')
            
            
            %% Part 2: For regression Analysis
            
            %% Downsampling for Brainstate analysis
            Theta_band = [6,12];Delta_band = [1.9,4];
            theta_freq = ((f>=Theta_band(1))&(f<Theta_band(2)));
            active_ts = t(:,TrackingVelocity>2);
            TrackingVelocity_active = TrackingVelocity(:,TrackingVelocity>2);
            active_theta  = (active_SPG(theta_freq,:));  % Isolate theta band from active SPG (for Pow/Vel & Frequ/Vel plots)
            %% Binning to 0.5sec
            ts_stepsize = t(3)-t(2); per_sec = 1/ts_stepsize;
            half_sec_elem = round(2*(numel(active_ts) / per_sec));
            if half_sec_elem ~= 0
                binned_t_halfsec = imresize(active_ts,[1,half_sec_elem],'bilinear');  %bins of 0.5sec
                binned_theta_halfsec = imresize(active_theta,[size(active_theta,1),half_sec_elem]); %bins of 0.5sec
                TrackingVelocity_active = imresize(TrackingVelocity_active,[1,half_sec_elem],'bilinear'); %bilinear, with bicubic Tracking velocity can become negative
            else
                binned_t_halfsec = [];  binned_theta_halfsec = []; TrackingVelocity_active = [];
            end
            %% Values
            [theta_amplit, theta_frequ_idx] = max(binned_theta_halfsec,[],1);
            theta_matrix = f(theta_freq);
            theta_max_frequ = theta_matrix(theta_frequ_idx);
            
            %% Half seconds binned data used for scatter plots (together s1, s2 and behaviour -> save hsb data of s1&s2 here and beh and make plot_Regression_pDay):
            hsb.TrackingVelocity = TrackingVelocity_active;
            hsb.theta_amplit = theta_amplit;
            hsb.theta_max_frequ = theta_max_frequ;
            
            save(fullfile(blockDir,sprintf('halfSecBinned_FreqPowVel.mat' )), 'hsb');
            clearvars hsb
            fprintf('Processing time for Session No.=%d: ',iii);
            toc(timer_session);
            
        end
        
    end
end
if ~isempty(failed_iii);     warndlg(sprintf('Failed days: %s', num2str(failed_iii)), 'Summary failed days');    end