function [] = runUnits()
clc
clearvars -except varargin
tic
addpath(genpath((pwd))); % need to modify this more specific
general_code_dir = fileparts(pwd);
addpath(fullfile(general_code_dir,'ImportantMetadata'));
addpath(fullfile(general_code_dir,'Tools'));
Silvia_PC_path
load(strcat(general_code_dir ,filesep, '+Figure8DataOrganization', filesep,'sessionInfo.mat'));

datadef2groups

for run_sessions = 1:3
    
    switch run_sessions
        case 1
            session = {'s1','s2'}; %  need to fix spikingproperties difference of intrinsic frequ to LFP (currently just using behavioural LFP)
        case 2
            session = {'d0','d2','d10'};  %or
        case 3
            session = {'d0','d2','d10','s1','s2'} %which is required for readECLR_silvia_2019.m
    end
    
    for age_gr = 1:2
        if age_gr == 1
            ageism = 'Young';
        elseif age_gr == 2
            ageism = 'Old';
        end
        for rec_area = 1:2
            if rec_area ==1
                recordinglayers = [1];
                name_reclayer = 'CA1';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
            elseif rec_area == 2
                recordinglayers = [2,3,4];
                name_reclayer = 'CA3DG';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
            end
            clearvars -except session ageism rec_area recordinglayers name_reclayer mainDirs sessInfo
            msgbox(sprintf('Starting processing runThetaAC for: %s - %s - %s', ageism, name_reclayer, strcat([session{:}])));
            
            mainDirs = {sessInfo.mainDir}';
            TFile = {sessInfo.tList}';
            analyzing_group = strcat(ageism, 'CA3rec' , name_reclayer);
            
            [group_A, group_B, animals_analysed, animal_no] = datadef2groups([1,2], recordinglayers, ageism);  %import which animals to analyse (search_geno, search_rec_lay, search_age)
            
            reanalyze = 0;
            
            plot_acfigures  = 0;
            
            
            all_corresp_anim = [];
            all_corresp_i = [];
            all_group = [];
            all_peak_fr_rate = [];
            all_avfrequ = [];
            all_burst_index = [];
            allWF_peakRatio = [];
            allWF_spikeWidth = [];
            all_intrins_frequ = [];
            all_thetaDiff2LFP = [];
            all_thetaMod_score = [];
            all_phase_mean = [];
            all_phase_pval = [] ;
            all_phase_std = [];
            all_ph_mlength = [];
            no_spikes_in_i= [];
            all_speedScore = [];
            all_speed_pVal = [];
            
            
            if reanalyze == 1;
                
                for group = 1:2
                    if group == 1, process_data = find(group_A);
                    elseif group == 2, process_data = find(group_B); end
                    
                    all_coll_ac = [];
                    
                    for i = process_data
                        fprintf('Processing, %s\n',num2str(i));

                        if ~exist(fullfile(mainDirs{i}, TFile{i}),'file')
                            errordlg(sprintf('TT file missing in i: %d %s',i,sessInfo(i).mainDir));
                            continue
                        end
                        
                        F = ReadFileList(fullfile(mainDirs{i}, TFile{i}));
                        if F{1} == -1  % jump if no cells
                            errordlg(sprintf('TT List empty in i: %d %s',i,sessInfo(i).mainDir));
                            continue
                        end
                        
                        %% Read spikes and define properties of all cells
                        
                        clearvars avg_frequ_cellsday peak_fr_rate interneurons
                        
                        [avg_frequ_cellsday, peak_fr_rate, max_burst, burst_index, intrins_frequ, theta_difference_to_LFP, thetaMod_score, coll_ac, ac_timing, phase_mean, phase_pval, phase_std, ph_mlength, speedScore, speed_pVal] = spiking_properties(sessInfo, i,F,session,plot_acfigures);
                        cellIDs = F;rec_day = i;
                        spk_properties_file = fullfile(sessInfo(i).mainDir, 'processedData', sprintf('spk_properties_%s.mat',cell2mat(session)));
                        save(spk_properties_file,'rec_day','cellIDs', 'avg_frequ_cellsday', 'peak_fr_rate', 'max_burst', 'burst_index', 'intrins_frequ', 'theta_difference_to_LFP',  'thetaMod_score', 'coll_ac', 'ac_timing', 'phase_mean', 'phase_pval', 'phase_std', 'ph_mlength','speedScore', 'speed_pVal')
                        

                        no_spikes_in_i = [no_spikes_in_i, i];
                        
                        %% Load waveformStats
                        wf_file = fullfile(sessInfo(i).mainDir, 'processedData', 'waveformStats.mat')
                        if exist(wf_file, 'file')
                            load(wf_file);
                            
                            if max([WF.PeakRatio])>5
                                for findweird_PRatio = 1:size(WF,2)
                                    if WF(findweird_PRatio).PeakRatio>5, errordlg(sprintf('Peak Ratio in %s %s: %s', num2str(i),WF(findweird_PRatio).cellID, WF(findweird_PRatio).PeakRatio)),end
                                end
                            end
                            if peak_fr_rate>200, errordlg(sprintf('Peak Ratio in %s %s: %s', num2str(i),WF(find(peak_fr_rate>200)).cellID, WF(find(peak_fr_rate>200)).PeakRatio)),end
                            
                            
                        else
                            continue
                        end
                        for allcc = 1:size(WF,2)
                            if ~(allcc == find(strcmp(WF(allcc).cellID,F))); %check to which cell in F the cellID corrsponds (typically they are both same sequence)
                                errordlg(sprintf('Mismatch in i: %s \n cell: %s', num2str(i), WF(allcc).cellID));
                                return
                            end
                        end
                        if exist('WF', 'var')  % only collect data if waveform of the cells are available
                            if size([WF.PeakRatio],2)~=size(avg_frequ_cellsday,1)
                                error('mismatch!!')
                            end
                            allWF_peakRatio = [allWF_peakRatio; [WF.PeakRatio]'];
                            allWF_spikeWidth = [allWF_spikeWidth; [WF.spikeWidth]'];
                            clear WF
                            
                            all_avfrequ = [all_avfrequ; avg_frequ_cellsday];
                            all_peak_fr_rate = [all_peak_fr_rate; peak_fr_rate];
                            all_burst_index = [all_burst_index; burst_index];
                            all_intrins_frequ = [all_intrins_frequ; intrins_frequ];
                            all_thetaDiff2LFP =  [all_thetaDiff2LFP; theta_difference_to_LFP];
                            all_thetaMod_score = [all_thetaMod_score; thetaMod_score];
                            all_coll_ac = cat(2,all_coll_ac,coll_ac);
                            
                            all_phase_mean = [all_phase_mean; phase_mean];
                            all_phase_pval = [all_phase_pval; phase_pval] ;
                            all_phase_std = [all_phase_std; phase_std];
                            all_ph_mlength = [all_ph_mlength; ph_mlength];
                            
                            all_speedScore = [all_speedScore; speedScore];
                            all_speed_pVal = [all_speed_pVal; speed_pVal];
                            
                            
                            all_corresp_anim(size(all_group,1)+1:size(all_group,1)+size(intrins_frequ,1),1) = sessInfo(i).animal; % use to check which animals
                            all_group(size(all_group,1)+1:size(all_group,1)+size(intrins_frequ,1),1) = group;
                            all_corresp_i(size(all_group,1)+1:size(all_group,1)+size(intrins_frequ,1),1) = i;
                            
                        end
                        %end
                        
                        
                        
                    end
                    
                    
                    
                    switch group
                        case 1
                            group1_ac = all_coll_ac;
                        case 2
                            group2_ac = all_coll_ac;
                    end
                    
                end

    all_coll_ac = [group1_ac, group2_ac];
    
    interneur = (all_avfrequ>5) .* (allWF_peakRatio<1.4);
    pyr = (all_avfrequ<10) .* (allWF_peakRatio>1.4);
    disp('Processing Time:');
    toc
    
    def1 = (all_avfrequ<10 & all_avfrequ>5 & allWF_peakRatio<1.1);
    def2 = (all_avfrequ<20 & all_avfrequ>10 & allWF_peakRatio<1.3);
    def3 = (all_avfrequ<30 & all_avfrequ>20 & allWF_peakRatio<1.5);
    def4 = ( all_avfrequ>30 & allWF_peakRatio<1.7);
    interneur2 = (def1 | def2 | def3 | def4);
    pyr2 = (all_avfrequ<10) .* (allWF_peakRatio>1.1);

                                                
                separation1 = (0.5 + (0.1* all_avfrequ)) - allWF_peakRatio;
                separation2 = (-1 + (0.4 * all_avfrequ)) - allWF_peakRatio; %interneurons >0 on separation line pyr<0
                pyr =  allWF_peakRatio>0.9 &  separation2<0 & all_avfrequ<10; %separation1<0
                
                if ismember(1,recordinglayers) %CA1 [1]
                    interneur =  ~pyr & all_avfrequ>5 & all_avfrequ<60 & allWF_peakRatio<2 & separation1>0;    
                else %CA3DG [2,3,4]
                    interneur =  ~pyr & all_avfrequ>5 & all_avfrequ<60 & allWF_peakRatio<2; 
                end
                
                
                fr = [0:60];
                
                peakR_graph = -1 + (0.4 * fr);
                figure
                plot(peakR_graph, fr,'k'), hold on
                if ismember(1,recordinglayers)  % only if it is CA1, the additional criterea is used
                    peakR_1 = 0.5 + (0.1* fr);
                    plot(peakR_1, fr,'k'),hold on
                end
                plot([0.9, 0.9],[0,10],'k')
                plot([2, 2],[10,60],'k')
                plot([0, 2],[5,5],'k')
                plot([0.9, 5],[10,10],'k')
                scatter(allWF_peakRatio(pyr&(all_group==1)),all_avfrequ(pyr&(all_group==1)),'b'),hold on
                scatter(allWF_peakRatio(pyr&(all_group==2)),all_avfrequ(pyr&(all_group==2)),'b*'),hold on
                scatter(allWF_peakRatio(interneur&(all_group==1)),all_avfrequ(interneur&(all_group==1)),'r')
                scatter(allWF_peakRatio(interneur&(all_group==2)),all_avfrequ(interneur&(all_group==2)),'r*')
                scatter(allWF_peakRatio(~interneur&~pyr&(all_group==1)),all_avfrequ(~interneur&~pyr&(all_group==1)),'k')
                scatter(allWF_peakRatio(~interneur&~pyr&(all_group==2)),all_avfrequ(~interneur&~pyr&(all_group==2)),'k*')
                xlim([0 5])
                ylim([0 60])
                hold off
                excluded = ~(pyr | interneur);
                unit_crit_folder = ('C:\Users\Silvia\Dropbox\Plots\UnitCriteria'); mkdir(unit_crit_folder);
                
                saveas(gcf, fullfile(unit_crit_folder,sprintf('UnitCriteria_%s_%s', analyzing_group, cell2mat(session))),'pdf');
                close all
                
                disp('Processing time, group');
                all_coll_ac = [group1_ac, group2_ac];

        X = cat(2,all_avfrequ, (allWF_peakRatio')); %scaling to equalize distance
        
        %%
        figure
        opts = statset('Display','final');
        [cidx, ctrs] = kmeans(X, 3, 'Distance','city', ...
            'Replicates',5, 'Options',opts);
        plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
        
        %% or
        c = clusterdata(X,'linkage','ward','savememory','on','maxclust',4)
        figure
        scatter(all_avfrequ(c==1),allWF_peakRatio(c==1),'r'),hold on
        scatter(all_avfrequ(c==2),allWF_peakRatio(c==2),'g'),hold on
        scatter(all_avfrequ(c==3),allWF_peakRatio(c==3),'b'),hold on
        scatter(all_avfrequ(c==4),allWF_peakRatio(c==4),'k'),hold on
                %}
                save(sprintf('AC_Data_11_12_2019_%s_%s.mat', analyzing_group, cell2mat(session)));
            else
                load(sprintf('AC_Data_11_12_2019_%s_%s.mat', analyzing_group, cell2mat(session)));
                %save(sprintf('AC_Data_08_02_2019_%s_%s.mat', analyzing_group, cell2mat(session)));
            end  % Finish reanalysis
            
                 
            %% Now export data per animal
            prismfolder = ('C:\Users\Silvia\Dropbox\SilviaProjectCode\PrismCSVFiles'); %mkdir(prismfolder);
            per_animalfolder = fullfile(prismfolder, 'Unit_stats_per_animal');
            mkdir(per_animalfolder);
            for subgrp = 1:2
                
            [pAnim_intr_frequ_pyr, pAnim_intr_frequ_int, pAnim_thetaMod_pyr,pAnim_thetaMod_int] = deal([]);
            [pAnim_thetaPhase_pyr, pAnim_thetaPhase_int, pAnim_thetaVectLength_pyr, pAnim_thetaVectLength_int, pAnim_speedscore_int, pAnim_BurstIdx_pyr, pAnim_BurstIdx_int] = deal([]);
            [pAnim_ave_firing_pyr, pAnim_ave_firing_int] = deal([]);
            per_animal_g = unique(all_corresp_anim(all_group==subgrp));
            subgrp_name = sprintf('Group%s',  num2str(subgrp));
            ticker_no=0;
            active_cells = all_avfrequ>0.1;
            for ticker_no = 1:numel(per_animal_g)
                anim_numb = per_animal_g(ticker_no);
                animals_corresp_indiv_numbers(ticker_no,1) = anim_numb; 
                %% Average Firing Rate
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group,'_', subgrp_name)), sprintf('Average Firing - Princ - Animal %s', num2str(anim_numb)), all_avfrequ((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Average Firing - Intern - Animal %s', num2str(anim_numb)), all_avfrequ((all_corresp_anim==anim_numb) & interneur)); % writing csv files

                pAnim_ave_firing_pyr(ticker_no,1) = nanmean(all_avfrequ((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_ave_firing_int(ticker_no,1) = nanmean(all_avfrequ((all_corresp_anim==anim_numb) & interneur));            

                %% Intrinsic Frequ
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Intrinsic Frequ - Princ - Animal %s', num2str(anim_numb)), all_intrins_frequ((all_corresp_anim==anim_numb) & pyr)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group,'_', subgrp_name)), sprintf('Intrinsic Frequ - Princ Active - Animal %s', num2str(anim_numb)), all_intrins_frequ((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Intrinsic Frequ - Intern - Animal %s', num2str(anim_numb)), all_intrins_frequ((all_corresp_anim==anim_numb) & interneur)); % writing csv files

                pAnim_intr_frequ_pyr(ticker_no,1) = nanmean(all_intrins_frequ((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_intr_frequ_int(ticker_no,1) = nanmean(all_intrins_frequ((all_corresp_anim==anim_numb) & interneur));
                %% Theta mod index
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Theta modulation score - Princ Active - Animal %s', num2str(anim_numb)), all_thetaMod_score((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Theta modulation score - Intern - Animal %s', num2str(anim_numb)), all_thetaMod_score((all_corresp_anim==anim_numb) & interneur)); % writing csv files
                pAnim_thetaMod_pyr(ticker_no,1) = nanmean(all_thetaMod_score((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_thetaMod_int(ticker_no,1) = nanmean(all_thetaMod_score((all_corresp_anim==anim_numb) & interneur));            

                %% Theta mean phase
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Theta mean phase - Princ Active - Animal %s', num2str(anim_numb)), all_phase_mean((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Theta mean phase - Intern - Animal %s', num2str(anim_numb)), all_phase_mean((all_corresp_anim==anim_numb) & interneur)); % writing csv files
                pAnim_thetaPhase_pyr(ticker_no,1) = nanmean(all_phase_mean((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_thetaPhase_int(ticker_no,1) = nanmean(all_phase_mean((all_corresp_anim==anim_numb) & interneur));            

                %% Theta mean vector length
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Theta mean phase - Princ Active - Animal %s', num2str(anim_numb)), all_ph_mlength((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Theta mean phase - Intern - Animal %s', num2str(anim_numb)), all_ph_mlength((all_corresp_anim==anim_numb) & interneur)); % writing csv files
                pAnim_thetaVectLength_pyr(ticker_no,1) = nanmean(all_ph_mlength((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_thetaVectLength_int(ticker_no,1) = nanmean(all_ph_mlength((all_corresp_anim==anim_numb) & interneur));            

                %% Speedscore  (only for interneurons)
                %write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Speedscore - Princ Active - Animal %s', num2str(anim_numb)), all_phase_mean((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Speedscore - Intern - Animal %s', num2str(anim_numb)), all_speedScore((all_corresp_anim==anim_numb) & interneur)); % writing csv files
                %pAnim_speedscore_pyr(ticker) = nanmean(all_speedScore((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_speedscore_int(ticker_no,1) = nanmean(all_speedScore((all_corresp_anim==anim_numb) & interneur));            

               %% Burst index
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Burst index - Princ Active - Animal %s', num2str(anim_numb)), all_burst_index((all_corresp_anim==anim_numb) & pyr & active_cells)); % writing csv files
                write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), sprintf('Burst index - Intern - Animal %s', num2str(anim_numb)), all_burst_index((all_corresp_anim==anim_numb) & interneur)); % writing csv files
                pAnim_BurstIdx_pyr(ticker_no,1) = nanmean(all_burst_index((all_corresp_anim==anim_numb) & pyr & active_cells));
                pAnim_BurstIdx_int(ticker_no,1) = nanmean(all_burst_index((all_corresp_anim==anim_numb) & interneur));            

            end
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Animal numbers', animals_corresp_indiv_numbers); % writing csv files
            
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Princip - Average Firing Rate', pAnim_ave_firing_pyr); % writing csv files
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur - Average Firing Rate', pAnim_ave_firing_int); % writing csv files

            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Princip - Intrinsic frequ', pAnim_intr_frequ_pyr); % writing csv files
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur - Intrinsic frequ', pAnim_intr_frequ_int); % writing csv files

            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Princip - Theta Mod Score', pAnim_thetaMod_pyr); 
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur - Theta Mod Score', pAnim_thetaMod_int); 

            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Princip - Theta Phase', pAnim_thetaPhase_pyr); 
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur - Theta Phase', pAnim_thetaPhase_int);             
                      
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Princip - Theta mean vect length', pAnim_thetaVectLength_pyr); 
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur - Theta mean vect lengt', pAnim_thetaVectLength_int);             
            
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur -Speedscore', pAnim_speedscore_int); 
           
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Princip - Burst index', pAnim_BurstIdx_pyr); 
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Average per animal - Interneur - Burst index', pAnim_BurstIdx_int);                           
            end  %end export for subgroup (WT vs TG of this condition)

        
        
            
            %% Write CSV Files
            prismfolder = ('C:\Users\Silvia\Dropbox\SilviaProjectCode\PrismCSVFiles'); mkdir(prismfolder);
            waveformfolder = fullfile(prismfolder, 'Thetacorrelation')
            mkdir(waveformfolder);
            
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Speed score',all_speedScore((all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Speed score',all_speedScore((all_group==2)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells  Speed score',all_speedScore((all_group==1)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells  Speed score',all_speedScore((all_group==2)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Theta modulation',all_thetaMod_score((all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Theta modulation',all_thetaMod_score((all_group==2)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells Theta modulation',all_thetaMod_score((all_group==1)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells Theta modulation',all_thetaMod_score((all_group==2)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Intrinsic frequency',all_intrins_frequ((all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Intrinsic frequency',all_intrins_frequ((all_group==2)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells Intrinsic frequency',all_intrins_frequ((all_group==1)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells Intrinsic frequency',all_intrins_frequ((all_group==2)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Intrinsic frequency -6Hz and 14Hz',all_intrins_frequ((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Intrinsic frequency -6Hz and 14Hz',all_intrins_frequ((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==2)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells Intrinsic frequency -6Hz and 14Hz',all_intrins_frequ((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==1)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells Intrinsic frequency -6Hz and 14Hz',all_intrins_frequ((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==2)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            
            
 
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Intrinsic freq Diff to LFP',all_thetaDiff2LFP((all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Intrinsic freq Diff to LFP',all_thetaDiff2LFP((all_group==2)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells Intrinsic freq Diff to LFP',all_thetaDiff2LFP((all_group==1)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells Intrinsic freq Diff to LFP',all_thetaDiff2LFP((all_group==2)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Intrinsic freq Diff to LFP -6Hz and 14Hz',all_thetaDiff2LFP((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Intrinsic freq Diff to LFP -6Hz and 14Hz',all_thetaDiff2LFP((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==2)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells Intrinsic freq Diff to LFP -6Hz and 14Hz',all_thetaDiff2LFP((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==1)&pyr&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells Intrinsic freq Diff to LFP -6Hz and 14Hz',all_thetaDiff2LFP((all_intrins_frequ>=6.5)&(all_intrins_frequ<=13.5)&(all_group==2)&pyr&(all_avfrequ>0.1)) ); % writing csv files

            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Aver firing rate notfilt',all_avfrequ((all_group==1)&interneur) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Aver firing rate notfilt',all_avfrequ((all_group==2)&interneur) ) ;
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Cells Aver firing rate notfilt',all_avfrequ((all_group==1)&pyr) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Cells Aver firing rate notfilt',all_avfrequ((all_group==2)&pyr) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Theta Phase Mean', all_phase_mean((all_group==1)&interneur&(all_avfrequ>0.1)) ); % writing csv files
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Theta Phase Mean', all_phase_mean((all_group==2)&interneur&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Theta Phase Mean', all_phase_mean((all_group==1)&pyr&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Theta Phase Mean', all_phase_mean((all_group==2)&pyr&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Theta Phase locking P-Value',  all_phase_pval((all_group==1)&interneur&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Theta Phase locking P-Value', all_phase_pval((all_group==2)&interneur&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Theta Phase locking P-Value', all_phase_pval((all_group==1)&pyr&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Theta Phase locking P-Value', all_phase_pval((all_group==2)&pyr&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Theta Phase vector length',  all_ph_mlength((all_group==1)&interneur&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Theta Phase vector length', all_ph_mlength((all_group==2)&interneur&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Theta Phase vector length', all_ph_mlength((all_group==1)&pyr&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Theta Phase vector length', all_ph_mlength((all_group==2)&pyr&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Theta Phase Standardev',  all_phase_std((all_group==1)&interneur&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Theta Phase Standardev', all_phase_std((all_group==2)&interneur&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Theta Phase Standardev', all_phase_std((all_group==1)&pyr&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Theta Phase Standardev', all_phase_std((all_group==2)&pyr&(all_avfrequ>0.1)) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Interneurons Burst Index',  all_burst_index((all_group==1)&interneur) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Interneurons Burst Index', all_burst_index((all_group==2)&interneur) );
            
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group1')), 'Principal Burst Index',  all_burst_index((all_group==1)&pyr&(all_avfrequ>0.1)) );
            write_csvfile(fullfile(waveformfolder, strcat(cell2mat(session),'_', analyzing_group, '_Group2')), 'Principal Burst Index', all_burst_index((all_group==2)&pyr&(all_avfrequ>0.1)) );
            

            
            
        end
    end
    
    
end

disp('Processing percentage of active cells');
%% Now loop again to calculate % of active cells where data is taken from d0d2d10 for number of active cells and d0d2d10s1s2 for number of all cells
for age_gr = 1:2
    if age_gr == 1
        ageism = 'Young';
    elseif age_gr == 2
        ageism = 'Old';
    end
    for rec_area = 1:2
        if rec_area ==1
            recordinglayers = [1];
            name_reclayer = 'CA1';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
        elseif rec_area == 2
            recordinglayers = [2,3,4];
            name_reclayer = 'CA3DG';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
        end
        clearvars -except session ageism rec_area recordinglayers name_reclayer mainDirs sessInfo
        analyzing_group = strcat(ageism, 'CA3rec' , name_reclayer); %YoungCA3recCA1 YoungCA3recCA3DG OldCA3recCA3DG
        
        prismfolder = ('C:\Users\Silvia\Dropbox\SilviaProjectCode\PrismCSVFiles'); %mkdir(prismfolder);
        per_animalfolder = fullfile(prismfolder, 'Unit_stats_per_animal');
        for subgrp = 1:2       
            load(sprintf('AC_Data_11_12_2019_%s_d0d2d10.mat', analyzing_group ));
            [anim_numb, per_animal_g, perc_active] = deal([]);
            per_animal_g = unique(all_corresp_anim(all_group==subgrp));
            subgrp_name = sprintf('Group%s',  num2str(subgrp));
            ticker_no=0; 
            for ticker_no = 1:numel(per_animal_g)
                anim_numb = per_animal_g(ticker_no);
                animals_corresp_indiv_numbers(ticker_no,1) = anim_numb;
                numb_active_cells(ticker_no,1) = sum((all_corresp_anim==anim_numb) & (all_avfrequ)>0.1 & pyr);
            end
            load(sprintf('AC_Data_11_12_2019_%s_d0d2d10s1s2.mat', analyzing_group ));
            for  ticker_no = 1:numel(per_animal_g)         
                anim_numb = per_animal_g(ticker_no);
                animals_corresp_indiv_numbers2(ticker_no,1) = anim_numb;
                numb_total_cells(ticker_no,1) = sum((all_corresp_anim==anim_numb) & pyr);               
            end
            perc_active =numb_active_cells ./ numb_total_cells;
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Active cells - Corresp Animal numbers', animals_corresp_indiv_numbers); % writing csv files
            write_csvfile(fullfile(per_animalfolder, strcat(cell2mat(session),'_', analyzing_group, '_', subgrp_name)), 'Active cells Percentage', perc_active); % writing csv files
        end
    end



end

close all
end

