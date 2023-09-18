function [avg_rate,max_fr_rate,allmax_sp_in_burst,burst_index,intrins_frequ, theta_difference_to_LFP, thetaMod_score,coll_ac,ac_timing, phase_mean, phase_pval, phase_std, ph_mlength, speedScore, speed_pVal] = spiking_properties(sessInfo, i,F,session,plot_acfigures)
% Function analyzes several spiking properties associated to cells:
% Avg Frequency, Peak firing rate, Max number of spikes in a burst, Burst
% index, Auto correlation / Theta modulation and a coresponding modulation
% score
addpath('C:\Users\Silvia\Dropbox\SilviaProjectCode');
all_phases = []; allph_corresp_ts = [];
ticker = 0;
for s = 1:numel(session)% sessInfo(i).sleepDirs
    ticker = ticker +1; %count current sessionnumber
    clear spikeTimes spikeData Pase eeg_ts_ds
    %% Read spikes
    spikeData = readSpikeDataOnly(fullfile(sessInfo(i).mainDir,session{s}) ,F);
    spikeTimes = fixSpikes(spikeData);
    collective_spikeTimes(:,s) = spikeTimes(:,1);
    
    %% Read theta phase file
    phasefile = fullfile(sessInfo(i).mainDir, session{s}, 'LFP', 'theta_phase.mat');
    if exist(phasefile, 'file') ~= 2
        thetaphase.run_getLFPphase(i, session(s)); % if the lfp theta phase file doesn't exist yet, create it now
    end
    load(phasefile, 'Phase', 'eeg_ts_ds');
    all_phases = [all_phases; Phase];  %used for phase locked
    allph_corresp_ts = [allph_corresp_ts; eeg_ts_ds];
    
    %% Read indata
    if session{s}(1) == 's'
        indata_s_file = fullfile(sessInfo(i).mainDir ,'processedData','indataS.mat');
        if exist(indata_s_file, 'file') ~= 2
            errordlg(sprintf('Error in indata sleep: %s', indata_s_file));
        end
        load(indata_s_file);
        s_num  = str2num(session{s}(2));
    elseif session{s}(1) == 'd'  %for behaviour data, look for d0, d2, d10
        load(fullfile(sessInfo(i).mainDir ,'processedData','indataB'));
        %% Get Theta Frequency here
        thetafrequ_file = fullfile(sessInfo(i).mainDir, 'processedData', 'theta_PSD_d0d2d10.mat');
          if exist(thetafrequ_file, 'file') ~= 2
         thetaphase.run_getPSDtheta(i, session(s)); % if the lfp theta phase file doesn't exist yet, create it now
         end
         load(thetafrequ_file, 'theta_frequ');
        
         switch session{s}(2)
            case '0'
                s_num  = 1;
            case '2'
                s_num  = 2;
            case '1'
                s_num  = 3;
        end
        
    else
        errordlg('Please define behav or sleep');
    end
    %total_time = tSp(end)- tSp(1);
    
    %% Readout: Average Firing rate
    %           Max no. of spikes in burst
    %           Mean burst index
    %           Max firing rate
    start_time = indata(s_num).t(1);
    end_time = indata(s_num).t(end);
    time_spent(ticker) = end_time-start_time;
    %ave_frequ = numel(tSp)/time_spent;
   
    
    for  ccc = 1:size(spikeTimes,1)
        clear tSp spikeintervall bursts start1 end1 burstlength spikes_per_burst
        tSp = spikeTimes{ccc};
   
        num_spikes(ccc,ticker) = numel(tSp); % keep track of number of spikes per cell per s1-d10        
        spikeintervall = diff(tSp);
        bursts = (spikeintervall < 0.006)';  % if spikeintervall is below 6 msec, it's a burst
        stst_indices = continuousRunsOfTrue(bursts);
        spikes_per_burst = 0;
        for bu = 1:size(stst_indices,1)
            inds = stst_indices(bu,1):stst_indices(bu,2)+1;  %add 1 to index since calculated from diff
            spikes_per_burst(bu) =  numel(inds);      
        end    
   
        max_sp_in_burst(ccc,ticker) = max(spikes_per_burst);  %collect for this cell and session max number of spikes per burst
        
       if nnz(bursts) == 0 
            max_sp_in_burst(ccc,ticker) = 0;   % if there are no bursts at all
            spikes_per_burst = 0;
        end
        
        burst_spikes(ccc,ticker) = sum(spikes_per_burst);     %collect for this cell and session total number of spikes that are part of bursts (for burst index)
        
        %%% Binning (using 0.25sec bins) for peak firing frequency
        for bins = 0.25:0.25:ceil(time_spent)
            timelimits= [(start_time+bins-0.25) (start_time+bins)];
            spikes_in_time(bins*4) = 4*nnz((spikeTimes{ccc}>timelimits(1)) .*(spikeTimes{ccc}<timelimits(2))); %revert the bins to 1 second steps so in Hz again
        end
        peak_fr_rate(ccc,ticker) = max(spikes_in_time);

        
        %% calculate speed modulation here
        %if strcmp(session{s}, 'd0')
        %disp('calculating speed score')
        %tic
        [speedScore_all(ccc,ticker), speed_pVal_all(ccc,ticker)] = getSpeedScore(indata(s_num).t,indata(s_num).v,tSp,2);
        %toc       
        %end
        
    end

    
end

spikeweighting = num_spikes ./ repmat(sum(num_spikes,2), [1,size(num_spikes,2)]);  % to adjust the speedScore by number of spikes in each behaviour/sleep session
speedScore = sum((spikeweighting .* speedScore_all), 2);
speed_pVal = sum((spikeweighting .* speed_pVal_all), 2);
%% Avg frequ
all_spikes = sum(num_spikes,2);
avg_rate = all_spikes/sum(time_spent);  %[in Hz]
%avg_rate_b = sum(num_spikes(:,1:3),2) / sum(time_spent(1,1:3));
%avg_rate_s = sum(num_spikes(:,4:5),2) / sum(time_spent(1,4:5));

%% Bursts
allmax_sp_in_burst = max(max_sp_in_burst,[],2);

% A burst index was defined as the ratio of spikes in bursts to all spikes.
summed_burst_spikes = sum(burst_spikes,2);
burst_index = summed_burst_spikes ./ all_spikes;

max_fr_rate = max(peak_fr_rate,[],2);


%% extract autocorrelation of all 5 sessions (d0-s2)
clear indata
coll_ac = []; coll_timing = [];
for  cc = 1:size(spikeTimes,1)   %autocorrelation of each cell
    tSp = [];
    for ses = 1:size(collective_spikeTimes, 2)
        tSp = [tSp; collective_spikeTimes{cc,ses}];
    end
    
    if numel(tSp) <= 10  %assign NaNs for cells that don't fire in the analyzed sessions
        intrins_frequ(cc,1) = NaN;
        thetaMod_score(cc,1) = NaN ;
        ac_nSp = NaN(121, 1);
        ac_timing = NaN(121, 1);
        
        phase_mean(cc,1) = NaN ;
        phase_pval(cc,1) = NaN ;
        phase_std(cc,1) = NaN ;
        ph_mlength(cc,1) = NaN ;
    else
        %% now run AC_Theta  (aka Intrinsic Frequency)
       % outdir = fullfile(sessInfo(i).mainDir, 'processedData', 'AutoCorrelations');  mkdir(outdir); img_outfile = strcat( F{cc}(1:end-2),'_',cell2mat(session));     
        % changed to look at autocorrel to select example
        outdir = ('C:\Users\Silvia\Desktop\CollectiveAutocorrelograms'); %mkdir(outdir);
        img_outfile = strcat('Animal_',num2str(sessInfo(i).animal),'_i_',num2str(i), '_', F{cc}(1:end-2),'_',cell2mat(session));             
        img_outfile = fullfile(outdir, img_outfile);
        
        [intrins_frequ(cc,1),thetaMod_score(cc,1), ac_timing, ac_nSp]  = AC_Theta(tSp',img_outfile, plot_acfigures);
        
        %% Phase analysis (see if Phase locked or not)
        phase_points = all_phases(nearestpoint(tSp, allph_corresp_ts));
        ph_po = circ_ang2rad(phase_points(~isnan(phase_points))); %get rid of NaNs and convert to rad
        phase_pval(cc,1) = circ_rtest(ph_po);
        phase_mean(cc,1) = circ_rad2ang(circ_mean(ph_po))+180;
        phase_std(cc,1) = circ_std(ph_po);
        ph_mlength(cc,1) = circ_r(ph_po);
        clear phase_points ph_po
    end
    coll_ac = cat(2, coll_ac, ac_nSp); %collect autocorrelograms of each cell
    
end
clear tSp collective_spikeTimes

theta_difference_to_LFP = intrins_frequ - theta_frequ;

if sum(abs([10, NaN, 3] - 5.5)>8)>=1
   error('Diff to high'); 
end

end

