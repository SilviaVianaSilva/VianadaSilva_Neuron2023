function [] = run_getLFPphase(varargin)
process_data = varargin{1};
session = varargin{2};
    
for iii = process_data
    clc
    fprintf(' -> Getting theta phase for i %s \n', num2str(iii));
    clearvars -except process_data iii failed_iii timer_prog session    
    load(fullfile('C:\Users\Silvia\Dropbox\SilviaProjectCode', '+Figure8DataOrganization','sessionInfo.mat'));
    for s = 1:numel(session)
        
        blockDir = fullfile(sessInfo(iii).mainDir, session{s});
        %% Load EEG
        channelInLayer = sprintf('CSC%d.ncs', sessInfo(iii).cellLayerChann);
        [eeg, eeg_fs,~] = readCRTsd(fullfile(blockDir, channelInLayer));
        eeg_raw = Data(eeg);%eeg_raw = eeg_raw(1:1000000,1);
        eeg_ts = Range(eeg) * 1e-4;
        
        eeg_fs_ds = eeg_fs / 16;
        eeg_raw_ds = downsample(eeg_raw, 16);
        eeg_ts_ds = downsample(eeg_ts, 16);
        [Phase, InstCycleFrequency, PerCycleFreq, signal_filtered] = thetaphase.DetectPhase(eeg_raw_ds,eeg_fs_ds);
  
        %% Save
        outfilename = fullfile(blockDir, 'LFP', 'theta_phase.mat');
        if exist(outfilename,'file'); delete(outfilename); end;
        if ~exist(fullfile(blockDir, 'LFP'),'dir'), mkdir(fullfile(blockDir, 'LFP')); end
        fprintf(' -> Saving Theta Phase of: i %s -- %s \n', num2str(iii), session{s});
        save(outfilename, 'Phase', 'InstCycleFrequency', 'PerCycleFreq', 'signal_filtered', 'eeg_ts_ds', 'eeg_fs_ds');
    
    end
end

end