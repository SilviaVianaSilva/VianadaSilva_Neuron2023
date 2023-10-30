% Determine Times per Area for LFP_Analysis in parsed maze
% And Extract regional EEGs
% Need to preprocessPathData before

disp(int2str(clock))
disp('Starting Extract regional EEGs')
clc

load(fullfile(codepath, '+Figure8DataOrganization', 'sessionInfo.mat'));
All_sessInfo = sessInfo; clear sessInfo

animal_numbers = unique([All_sessInfo.animal]);
failed_iii = [];

timer_prog = clock;
problematicSessions = {};
for iii = process_data
    if ~isempty(All_sessInfo(iii).animal)
        
        fprintf('\n -- Starting to process i number: %d --\n', iii);
        close all
        clearvars -except iii process_data All_sessInfo sessInfo error_in_i timer_prog aaa animal_numbers failed_iii problematicSessions  
        sessInfo = All_sessInfo(iii);
        timer_anim = tic;

        all_rois = {'return';'delay';'stem';'choice';'reward'};
        
        all_groups = {'WtCA3';'TgCA3';'WtCA1';'TgCA1';'WtDG';'TgDG'}; %replaced switch module
        sessInfo.group = all_groups{sessInfo.genotype};
        sessGroups = {sessInfo.group};
        groups = unique(sessGroups);
        blocks = {sessInfo.sessDirs};
       [~, linBnds] = fig8trialtemplate();
        ptempl = parsingtemplate('fig8mouse:rdscr');
        allBlockNames = unique([sessInfo.sessDirs]);
        avgRates = [];
        zoneOccupancy = [];
        g = 1;  	i = 1;
        group = groups{g};

        try
            [trialInfo, parsingInfo, pathData, pathDataIdeal, pathDataLin] = plotTrialPerCell.loadInfo(sessInfo, i);
        catch
            problematicSessions = [problematicSessions; sessInfo(i).mainDir];
            errordlg(sprintf('Error in linearizing path of %s', num2str(iii)));
            continue;
        end
   
        display(['About to do session ' sessInfo(i).mainDir]);
        for block = sessInfo(i).sessDirs
            disp('Analyzing: ');disp(block);
            clearvars spect eeg channelInLayer
            blockDir = fullfile(sessInfo(i).mainDir, block{1});
            LFP_dir = fullfile(blockDir,'LFP'); if ~exist(LFP_dir,'dir'), mkdir(LFP_dir),end
            tinfo = trialInfo.(block{1});
            pinfo = parsingInfo.(block{1});
            pdata = pathData.(block{1});
            pideal = pathDataIdeal.(block{1});
            plin = pathDataLin.(block{1});    
            tridx = trialInfo2simple.trial_startend_ind(tinfo);
            tridx = tridx(~tinfo.degen, :);
            trtype = plotTrialPerCell.categtable(tinfo);
              
            %% Start processing EEGs
            channelInLayer = sprintf('CSC%d.ncs', sessInfo.cellLayerChann); % Picks the channel that is in the layer
            eeg = readCRTsd(fullfile(blockDir, channelInLayer));
            frequs = [2:300] ;
            clearvars velocity successful all_regions_eegx all_regions_eegxidx all_regions_times all_regions_tdomainbin
            clearvars spect frequ ROI_duration_time vid_ts
            
            opt.wavelet.dsrate = 16;  %Used for spectogramm
            opt.wavelet.smoothingWin = 20;opt.wavelet.freqRes = 20;
            fprintf('\nRunning wavelet ');
            for trial = 1:length(tridx)
                for region = 1:5
                    roixy = ptempl(strcmp({ptempl.zone}, all_rois{region}));
                    trialx = [tridx(trial,1) tridx(trial,2)];
                    [eegx_region, eegidx_region, eegtime_region, tdomain_bin] = plotTrialPerCell.extractEEGEpoch(pideal, trialx, roixy, eeg, [0, 0]);  %isolate eeg from this area
                    if ~isstruct(eegx_region)  % animal not tracked in this part of the maze
                        errordlg(sprintf('Could extract EEG in region %d, lap %d animal/i %d, will erase this lap', region, trial, iii));
                        ROI_duration_time(trial,region) = NaN;
                        continue
                    else % regular case
                        times = eegtime_region(2) - eegtime_region(1);
                        all_regions_eegx{trial,region} = eegx_region;
                        all_regions_eegxidx{trial,region} = eegidx_region;
                        all_regions_times{trial,region} = times;
                        all_regions_tdomainbin{trial,region} = tdomain_bin;
                        velocity{trial,region} = pdata.v(tdomain_bin.startidx:tdomain_bin.endidx);
                        vid_ts{trial,region} = pdata.t(tdomain_bin.startidx:tdomain_bin.endidx);
                        fprintf('.');
                        [SPG, t, f, bandSpecgramFun] = specgramwwd(eegx_region.d,eegx_region.Fs, 2, 300,opt.wavelet);
                        spect{trial,region} = SPG;
                        frequ{trial,region} = f;
                        ROI_duration_time(trial,region) = eegtime_region(2) - eegtime_region(1);
                        
                        
                    end

                end % For this ROI
            end
            fprintf('.\n\n');
            %% Now delete trials that have missing data
            delete_row = isnan(mean(ROI_duration_time,2));
            all_regions_eegx(find(delete_row),:) = [];
            all_regions_eegxidx(find(delete_row),:) = [];
            all_regions_times(find(delete_row),:) = [];
            all_regions_tdomainbin(find(delete_row),:) = [];
            velocity(find(delete_row),:) = [];
            vid_ts(find(delete_row),:) = [];
            spect(find(delete_row),:) = [];
            frequ(find(delete_row),:) = [];
            
            successful = tinfo.success(~delete_row);%(1:size(ROI_duration_time,1));
            ROI_duration_time(find(delete_row),:) = [];
            
            %% Save
            filename = fullfile(LFP_dir,sprintf('%s_Animal_i%03d_%s.mat', group, iii,block{1}));
            save(filename,'sessInfo','ROI_duration_time','successful','spect','frequ','all_regions_eegxidx','all_regions_times','all_regions_tdomainbin','velocity', 'vid_ts');
            %% Save info in metafile for analysis
            metafilename=('.\metafile_RegionalEEGs.mat');
            if exist(metafilename,'file'), load(metafilename);end
            datainfo.group = group;
            datainfo.animal = sessInfo.animal;datainfo.age = sessInfo.age; datainfo.session = sessInfo.session;     datainfo.filename = filename;
            metainfo.(block{1})(iii) = datainfo;        
            save(metafilename,'metainfo');
        end % for block
        fprintf('Processing time for Animal No.=%d: ',iii);
        toc(timer_anim);
    end 
    
end  % end of cycling through days
elapsed_min = round(etime(clock,timer_prog)/60);
msgbox(sprintf('Total runtime: %d min\n',elapsed_min),'Runtime');