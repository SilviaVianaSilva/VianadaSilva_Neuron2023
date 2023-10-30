% Reads saved files from ConvertRegionalEEGs
%
%---------------------------------------------------------
% Matthias Haberl, 2018
%---------------------------------------------------------

clear

%%
for age_gr = 1:2
    if age_gr == 1
        ageism = 'young';
    elseif age_gr == 2
        ageism = 'old';
    end
    
    
    for rec_area = 1:2
        if rec_area ==1
            recordinglayers = [1];
            name_reclayer = 'CA1';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
        elseif rec_area == 2
            recordinglayers = [2,3,4];
            name_reclayer = 'DGCA3';  % 'DGCA3' 'CA1'  define here which name to use for your plots, can be different from all the layers accepted for recording
        end
        
        clearvars -except age_gr ageism rec_area recordinglayers name_reclayer
        accepted_layers = recordinglayers;
        fprintf('\n\n------------------->>> Processing: %s  ---  %s  <<<-------------------\n\n', ageism, name_reclayer);
        
        addpath('.\Metadata\');
        analyze_groups = 1:2;
        name_exprlayer = 'CA3';
        
        array_layers = {'CA1','CA3','DG','DG_CA3'};
        
        [group_A, group_B, animals_analysed, animal_no] = datadef2groups([1,2], accepted_layers, ageism);
        animals_groupA =  unique(animal_no(group_A));
        animals_groupB =  unique(animal_no(group_B));
        animals = {animals_groupA; animals_groupB};
        
        
        %% Load updated binning per time
        frequ_analysis.load_TimeBins
        x_length = sum(opt.binninByTime.d10);
        
        
        block = {'d0','d2','d10'};
        infolder = ('.\LFP_Analysis\AnimalSpects');
        
        
        plots_dir = fullfile(strcat('.\Plots\SleepNorm_',ageism, '_expr', name_exprlayer, '_rec' , name_reclayer,'\'));mkdir(plots_dir);
        plots_dir_animals = fullfile(plots_dir, 'perAnimal');mkdir(plots_dir_animals);
        analysis_dir = plots_dir;mkdir(analysis_dir); % Excel File Output
        
        csvtopdir = '.\CSVFiles\LFP_Analysis'; mkdir(csvtopdir);
        csvdir = fullfile(csvtopdir, 'sleepNormLFP'); mkdir(csvdir);
        
        %% Settings
        doplot_individual_animals = 0;
        doplot_powerFrequ_regions = 0;
        options.y_scale = [0.5 2];  options.colorcode=('jet');
        options.yaxis=[2 300];  options.xaxis = [0 x_length];
        
        for ddd = 1:3
            csv_delaydir = fullfile(csvdir ,block{ddd}); mkdir(csv_delaydir);
            for gg  = analyze_groups %group                
                bins = opt.binninByTime.(block{ddd});
                bounds = frequ_analysis.bins2bounds(bins);
                options.bounds= bounds;
                aver_spect = [];
                group_spect = [];
                group_spect_fails = [];
                combined_speed = [];
                combined_speed_fails = [];
                animals_quiet_power_spect = [];
                animals_rawpower_spect = [];
                animal_count = 0;
                [all_anim, all_delt_fr, all_theta_fr, all_beta_fr, all_gamma_fr, all_higamma_fr] = deal([]); %accum data for the csvfiles
                [all_delt_pow, all_theta_pow, all_beta_pow, all_gamma_pow, all_higamma_pow] = deal([]); %accum data for the csvfiles
                for aaa = animals{gg}
                    fprintf('==> Analyzing Group %s <==\n', num2str(gg));
                    animal_count = animal_count + 1;
                    tic
                    fprintf('-> Analysing: %s --- Animal: %s\n', block{ddd}, num2str(aaa));
                    %filename = fullfile(infolder,sprintf('Animal_%d_adjusted_spect_behav_%s.mat',aaa,block{ddd}));
                    
                    filename2 = [];
                    for LLL = 1:numel(accepted_layers)  % if more than 1 layer is accepted, circulate to find which name/layer was entered in Fig8 excel sheet
                        filename = fullfile(infolder,sprintf('Animal_%d_adjusted_spect_behav_recLay_%s_%s.mat',aaa,array_layers{accepted_layers(LLL)},block{ddd}));
                        if exist(filename,'file')
                            filename2 = filename;
                        end
                    end
                    if isempty(filename2)
                        errordlg(sprintf('No file found for: Animal_%d_adjusted_spect_behav_recLay_%s_%s.mat\n%s',aaa,array_layers{accepted_layers(LLL)},block{ddd}), infolder);
                    end
                    
                    if exist(filename2,'file')
                        load(filename2,'animal_spect','animal_sleepNormalized','f','animal_fails_sleepNormalized','animal_speed', 'animal_speed_fail','animal_quiet_power_spect');
                        
                        animals_rawpower_spect = cat(3,animals_rawpower_spect,mean(mean(animal_spect, 2),3));
                        animals_quiet_power_spect = cat(3,animals_quiet_power_spect,mean(animal_quiet_power_spect, 3)); % to plot power over freq
                        
                        group_spect = cat(3,group_spect,nanmean(animal_sleepNormalized,3));
                        group_spect_fails = cat(3,group_spect_fails,nanmean(animal_fails_sleepNormalized,3));
                        combined_speed = cat(1,combined_speed,nanmean(animal_speed,3));
                        combined_speed_fails = cat(1,combined_speed_fails,nanmean(animal_speed_fail,3));
                        
                        power_freq_animal = nanmean(animal_sleepNormalized,2);
                        
                        coll_frequ = NaN(size(power_freq_animal, 3), 5, 5);
                        coll_power = NaN(size(power_freq_animal, 3), 5, 5);
                        for rd = 1:size(power_freq_animal, 3)  %circulate through recording days
                            
                            %% Parsed
                            export_areas_to_excel = [1:6];
                            area_names = {'return', 'delay', 'stem', 'choice', 'reward','centr_arm'}; % 6th region stem+choice = central arm
                            excelfilename = fullfile(analysis_dir,sprintf('AllGroups_%s_%s.xls',block{ddd},date));
                            
                            for this_area = export_areas_to_excel
                                if this_area == 6
                                power_freq_animal_area =  mean(animal_sleepNormalized(:,bounds(3)+1:bounds(5),rd),2); % 6th region stem+choice = central arm    
                                else
                                power_freq_animal_area =  mean(animal_sleepNormalized(:,bounds(this_area)+1:bounds(this_area+1),rd),2);
                                end
                                [powers , indices, W,P] = findpeaks(power_freq_animal_area,'MinPeakProminence',0.01);
                                frequpks = f(indices);
                                
                                band_names = {'Delta' , 'Theta' , 'Beta', 'Gamma', 'HGamma'};
                                f_bands = {[2,6], [6, 12], [14, 26], [26, 50], [50, 150]};
                                
                                for bndno = 1:size(f_bands,2)
                                    f_range = f_bands{bndno};
                                    the_pks = find((frequpks>f_range(1)).*(frequpks<f_range(2)));
                                    if ~isempty(the_pks)
                                        [~,pk_idx] = max(powers(the_pks));
                                        coll_frequ(rd,this_area,bndno) = frequpks(the_pks(pk_idx));
                                        coll_power(rd,this_area,bndno) = powers(the_pks(pk_idx));
                                    else
                                        coll_frequ(rd,this_area,bndno) = NaN;
                                        coll_power(rd,this_area,bndno) = NaN;
                                    end
                                    
                                end
                                
                            end
                            
                        end  %finish circling through days
                        
                        for this_area = export_areas_to_excel
                            d = {'Animal Number', 'Delta-Multiple of Sleep', 'Theta-Multiple of Sleep', 'Beta-Multiple of Sleep ','Gamma Multiple','High-Gamma Multiple','Beta-Frequency (2-6)','Theta-Frequ (6-12)', 'Beta-Frequ (14-30)','Gamma-Frequ (30-60)','HighGamma-Frequ (60-120)'};
                            if animal_count == 1
                                xlswrite(excelfilename, d, sprintf('Group %s - Area %s', num2str(gg), area_names{this_area}), 'A1');
                            end
                            data2 = [aaa, nanmean(coll_power(:,this_area,1)), nanmean(coll_power(:,this_area,2)), nanmean(coll_power(:,this_area,3)),nanmean(coll_power(:,this_area,4)), nanmean(coll_power(:,this_area,5)), nanmean(coll_frequ(:,this_area,1)), nanmean(coll_frequ(:,this_area,2)), nanmean(coll_frequ(:,this_area,3)), nanmean(coll_frequ(:,this_area,4)),nanmean(coll_frequ(:,this_area,5))];
                            cellorganiz= {'A','B','C','D','E','F','G','H','I','J','K','L','M','O','P','Q'};
                            if size(data2,2)==11
                                xlswrite(excelfilename, data2, sprintf('Group %s - Area %s', num2str(gg), area_names{this_area}), [cellorganiz{1},num2str(animal_count+1)]);
                            else
                                errordlg('Error finding peaks in animal %s', num2str(aaa))
                            end
                            clear data2
                        end
                        
                        %explanation: coll_frequ(recording days, parsed area, frequ band)
                        all_anim = cat(1, all_anim, aaa);
                        all_delt_fr = cat(1, all_delt_fr, nanmean(coll_frequ(:,:,1),1));
                        all_theta_fr = cat(1, all_theta_fr, nanmean(coll_frequ(:,:,2),1));
                        all_beta_fr = cat(1, all_beta_fr, nanmean(coll_frequ(:,:,3),1));
                        all_gamma_fr = cat(1, all_gamma_fr, nanmean(coll_frequ(:,:,4),1));
                        all_higamma_fr = cat(1, all_higamma_fr, nanmean(coll_frequ(:,:,5),1));
                        all_delt_pow  = cat(1, all_delt_pow, nanmean(coll_power(:,:,1),1));
                        all_theta_pow = cat(1, all_theta_pow, nanmean(coll_power(:,:,2),1));
                        all_beta_pow = cat(1, all_beta_pow, nanmean(coll_power(:,:,3),1));
                        all_gamma_pow = cat(1, all_gamma_pow, nanmean(coll_power(:,:,4),1));
                        all_higamma_pow = cat(1, all_higamma_pow, nanmean(coll_power(:,:,5),1));
                        clearvars coll_frequ coll_power
                        
                        
                        %cellorganiz= {'A','B','C','D','E','F','G','H','I','J','K','L'};
                        %for ccc = 1:12  % takes care that data are written always in the right cells, even if there are empty cells
                        %    if ~isempty(data2{ccc})
                        %        xlswrite(excelfilename, data2{ccc}, sprintf('Group %s', num2str(gg)), [cellorganiz{ccc},num2str(animal_count)]);
                        %    end
                        %end
                        %xlswrite(excelfilename, low_gamma_frequ_1, num2str(aaa), 'B2');
                        %xlswrite(excelfilename, gamma_frequ_1, num2str(aaa), 'C2');
                        %xlswrite(excelfilename, Theta_frequ_2, num2str(aaa), 'D2');
                        %xlswrite(excelfilename, low_gamma_frequ_2, num2str(aaa), 'E2');
                        %xlswrite(excelfilename, gamma_frequ_2, num2str(aaa), 'F2');
                        % xlswrite(excelfilename, PKS', num2str(aaa), 'B1');
                        % xlswrite(excelfilename, LOCS, num2str(aaa), 'B2');
                        % xlswrite(excelfilename, W, num2str(aaa), 'B3');
                        % xlswrite(excelfilename, P', num2str(aaa), 'B4');
                        
                        
                        %% per Animal image
                        if doplot_individual_animals == 1
                            figure_title = sprintf('Group %d -- SleepNormalized -- %s -- Animal %d', gg, block{ddd}, aaa);
                            filename = fullfile(plots_dir_animals,sprintf('Aver_sleepNorm_spect_Group_%d_%s_Animal_%d_log', gg, block{ddd},aaa));
                            frequ_analysis.plot_spectogr(nanmean(animal_sleepNormalized,3),f,figure_title,filename,options);
                            
                            
                            figure_title = sprintf('Group %d -- SleepNormalized -- %s -- Animal %d -- failed trials', gg, block{ddd}, aaa);
                            filename = fullfile(plots_dir_animals,sprintf('Aver_sleepNorm_spect_Group_%d_%s_Animal_%d_log_failedtrials', gg, block{ddd},aaa));
                            frequ_analysis.plot_spectogr(nanmean(animal_fails_sleepNormalized,3),f,figure_title,filename,options);
                        end
                        
                    end
                    toc
                    clearvars animal_spect animal_sleepNormalized animal_fails_sleepNormalized animal_speed animal_speed_fail animal_quiet_power_spect
                end
                %% export into csvfile here
                genotype = {'wt', 'tg'};
                csv_subdir = fullfile(csv_delaydir, sprintf('%s_%s_expr%s_rec%s', genotype{gg}, ageism, name_exprlayer, name_reclayer)); mkdir(csv_subdir);
                write_csvfile(csv_subdir, sprintf('Animal delta frequ p region'), all_delt_fr');
                write_csvfile(csv_subdir, sprintf('Animal theta frequ p region'), all_theta_fr');
                write_csvfile(csv_subdir, sprintf('Animal beta frequ p region'), all_beta_fr');
                write_csvfile(csv_subdir, sprintf('Animal gamma frequ p region'), all_gamma_fr');
                write_csvfile(csv_subdir, sprintf('Animal high gamma frequ p region'), all_higamma_fr');
                write_csvfile(csv_subdir, sprintf('Animal delta multiple of sleep p region'), all_delt_pow');
                write_csvfile(csv_subdir, sprintf('Animal theta multiple of sleep p region'), all_theta_pow');
                write_csvfile(csv_subdir, sprintf('Animal beta multiple of sleep p region'), all_beta_pow');
                write_csvfile(csv_subdir, sprintf('Animal gamma multiple of sleep p region'), all_gamma_pow');
                write_csvfile(csv_subdir, sprintf('Animal high gamma multiple of sleep p region'), all_higamma_pow');
                write_csvfile(csv_subdir, 'AnimalNumbers', all_anim');
                
                   %% Now export individual frequ per area into different csv files
                for this_area = export_areas_to_excel
                    csv_subdir = fullfile(csv_delaydir, sprintf('%s_%s_expr%s_rec%s', genotype{gg}, ageism, name_exprlayer, name_reclayer),sprintf('Area_%s', area_names{this_area}) ); mkdir(csv_subdir);
                    write_csvfile(csv_subdir, sprintf('Animal delta frequ %s', area_names{this_area}), all_delt_fr(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal theta frequ  %s', area_names{this_area}), all_theta_fr(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal beta frequ  %s', area_names{this_area}), all_beta_fr(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal gamma frequ %s', area_names{this_area}), all_gamma_fr(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal high gamma frequ %s', area_names{this_area}), all_higamma_fr(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal delta power %s', area_names{this_area}), all_delt_pow(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal theta power %s', area_names{this_area}), all_theta_pow(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal beta power %s', area_names{this_area}), all_beta_pow(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal gamma power %s', area_names{this_area}), all_gamma_pow(:,this_area)');
                    write_csvfile(csv_subdir, sprintf('Animal high gamma power %s', area_names{this_area}), all_higamma_pow(:,this_area)');
                    write_csvfile(csv_subdir, 'AnimalNumbers', all_anim');
                end

                clearvars all_anim all_delt_fr all_theta_fr all_beta_fr all_gamma_fr all_higamma_fr
                
                aver_spect = nanmean(group_spect,3);
                aver_fail_spect  = nanmean(group_spect_fails,3);
                outfilename = fullfile(infolder,sprintf('Group_sleepNorm_%d_%s', gg, block{ddd}));
                save(outfilename,'aver_spect','aver_fail_spect','f');
                
                %% To plot power/frequ of beh and sleep1 of WT and TG
                power_frequ_beh{gg} = animals_rawpower_spect;
                power_frequ_sleep{gg} = animals_quiet_power_spect;
                
                %% Isolate power/frequency plot per group, normalized to sleep
                if doplot_powerFrequ_regions
                    % return: red
                    % del: black
                    % stem: yellow
                    % choice: blue
                    % reward: purple
                    % centr arm: green                  
                    line_color=  {[1 0 0];[0.3144 0.314 0.314];[1 0.60 0]; ...
                        [0 0.447 0.741]; [0.5451 0 0.5451];[0 0.6 0]}; %centr arm: green reward now: orange
                    legendlabel= {'return','delay','stem','choice','reward','centr_arm'};
                    for area = 1:6
                        if area == 6
                        power_freq_group_onearea =  mean(aver_spect(:,bounds(3)+1:bounds(5)),2);
                        else
                        power_freq_group_onearea =  mean(aver_spect(:,bounds(area)+1:bounds(area+1)),2);
                        end
                        parsed{area} = plot(f, power_freq_group_onearea,'LineWidth',1 ,'Color',(line_color{area}));hold on
                        title(sprintf('Group %d -- SleepNormalized -- %s -- All', gg, block{ddd}));
                        %savefig(filename)
                        set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                        xlabel('Freq (Hz)'); ylabel('Power normalized to sleep');
                        %print('-dpdf', fullfile(plots_dir,sprintf('Aver_sleepNorm_spect_Group_%d_%s_log_parsed_%s', gg, block{ddd},legendlabel{area})))
                        
                    end
                    hold off;
                    print('-dpdf', fullfile(plots_dir,sprintf('Aver_sleepNorm_spect_Group_%d_%s_log_parsed_ALL', gg, block{ddd})));
                    print('-dpng', fullfile(plots_dir,sprintf('Aver_sleepNorm_spect_Group_%d_%s_log_parsed_ALL', gg, block{ddd})));
                end
                %% Compare power/frequ between correct and incorrect trials
                clear power_freq_group_onearea
                figure
                line_color=  {[0 0 0];[1 0 0];[0.1 0.1 0.1];[1 0.2 0.2];[1 0 0];[0 0.600000023841858 0.6];};
                power_frequ_group_correct = nanmean(group_spect,2);
                power_frequ_group_incorrect = nanmean(group_spect_fails,2);
                plot(f, mean(power_frequ_group_correct,3),'LineWidth',1 ,'Color',(line_color{gg}));hold on
                plot(f, mean(power_frequ_group_incorrect,3),'--','LineWidth',1 ,'Color',(line_color{gg}));hold on
                title(sprintf('Group %d -- Correct vs incorrect -- %s -- All Regions', gg, block{ddd}));
                %savefig(filename)
                set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                xlabel('Freq (Hz)'); ylabel('Power in behav divided by sleep');ylim(gca,[0.2 2]);
                print('-dpdf', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_Group_%d_%s_log_allRegions', gg, block{ddd})))
                print('-dpng', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_Group_%d_%s_log_allRegions', gg, block{ddd})))
                close all
                
                line_color=  {[1 0 0];[0.314 0.314 0.314];[1 0.6 0]; ...
                     [0 0.447 0.741];[0.5451 0 0.5451];[0 0.6 0]}; %centr arm: green reward now: orange};
                legendlabel= {'return','delay','stem','choice','reward','centr_arm'};
                for area = 1:6
                    figure
                    if area == 6
                    power_freq_group_area{area,gg} =  nanmean(group_spect(:,bounds(3)+1:bounds(5),:),2);
                    power_freq_group_area_incorr{area,gg} =  nanmean(group_spect_fails(:,bounds(3)+1:bounds(5),:),2);
                    else
                    power_freq_group_area{area,gg} =  nanmean(group_spect(:,bounds(area)+1:bounds(area+1),:),2);
                    power_freq_group_area_incorr{area,gg} =  nanmean(group_spect_fails(:,bounds(area)+1:bounds(area+1),:),2);
                    end
                    plot(f, mean(power_freq_group_area{area,gg},3),'LineWidth',1 ,'Color',(line_color{area}));hold on
                    plot(f, mean(power_freq_group_area_incorr{area,gg},3),'--','LineWidth',1 ,'Color',(line_color{area}));
                    title(sprintf('Group %d -- SleepNormalized -- %s -- All', gg, block{ddd}));hold off
                    %savefig(filename)
                    set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                    xlabel('Freq (Hz)'); ylabel('Power normalized to sleep');ylim(gca,[0.2 2]);
                    print('-dpdf', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_Group_%d_%s_log_parsed_%s', gg, block{ddd},legendlabel{area})))
                    print('-dpng', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_Group_%d_%s_log_parsed_%s', gg, block{ddd},legendlabel{area})))
                    close
                end
                                
                f_linear = sort(min(f):0.1:max(f),'descend');
                spect_linear = interp1(f,aver_spect,f_linear,'nearest');
                
                %% PLot Spectograms here
                
                figure_title = sprintf('Group %d -- SleepNormalized -- %s', gg, block{ddd});
                
                filename = fullfile(plots_dir,sprintf('Aver_sleepNorm_spect_Group_%d_%s_log', gg, block{ddd}));
                frequ_analysis.plot_spectogr_MH(aver_spect,f,figure_title,filename,options); % Turn off to be quicker
                
                figure_title = sprintf('Group %d -- SleepNormalized -- %s -- failed trials', gg, block{ddd});                
                filename = fullfile(plots_dir,sprintf('Aver_sleepNorm_spect_Group_%d_%s_log_failed_trials', gg, block{ddd}));
                frequ_analysis.plot_spectogr_MH(aver_fail_spect,f,figure_title,filename,options); % Turn off to be quicker
                
                %% Pooling velocity into speed per group
                group_speed{gg} = combined_speed;
                group_speed_fails{gg} = combined_speed_fails;
                
            end
            
            for area = 1:5
                line_color=  {[0 0 0];[1 0 0];[0 0 0];[0.2 1 0];[0 0 0];[0 0.2 0.9]};
                legendlabel= {'return','delay','stem','choice','reward'};
                figure
                for gg = analyze_groups
                    plot(f, mean(power_freq_group_area{area,gg},3),'LineWidth',1 ,'Color',(line_color{gg}));hold on
                    plot(f, mean(power_freq_group_area_incorr{area,gg},3),'--','LineWidth',1 ,'Color',(line_color{gg}));hold on
                end
                title(sprintf('Area %d -- SleepNormalized -- %s -- All', area, block{ddd}));hold off
                %savefig(filename)
                set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                xlabel('Freq (Hz)'); ylabel('Power normalized to sleep');ylim(gca,[0.2 2.5]);
                print('-dpdf', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_AllGenotype_%s_log_parsed_%s', block{ddd},legendlabel{area})))
                print('-dpng', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_AllGenotype_%s_log_parsed_%s', block{ddd},legendlabel{area})))
            end
            
            %% Plot speed for WT/TG old animals (group 1 and 2)
            addpath('Tools');
            genotypes = {'wt','tg'};
            figure
            graph_WT_correct = shadedErrorBar([], nanmean(group_speed{analyze_groups(1)},1),nanstd(group_speed{analyze_groups(1)}), 'k-',1);hold on
            ylim(gca,[0 30]);
            set(gca,'XTick',bounds); %,'XTickLabel',{'return','delay','stem','choice','reward'});
            print('-dpdf', fullfile(plots_dir,sprintf('Velocity_%s_%s_%s_%s',ageism, name_reclayer, block{ddd}, genotypes{1})),'-painters');
 
            close
            figure
            graph_TG_correct = shadedErrorBar([], nanmean(group_speed{analyze_groups(2)},1),nanstd(group_speed{analyze_groups(2)}), 'g-',1);hold off            
            ylim(gca,[0 30]);
            set(gca,'XTick',bounds); %,'XTickLabel',{'return','delay','stem','choice','reward'});
            print('-dpdf', fullfile(plots_dir,sprintf('Velocity_%s_%s_%s_%s',ageism, name_reclayer, block{ddd}, genotypes{2})),'-painters');
  
            close
            figure
            graph_WT_incorrect = shadedErrorBar([], nanmean(group_speed_fails{analyze_groups(1)}),frequ_analysis.sem(group_speed_fails{analyze_groups(1)}), 'k-',1);hold on
            ylim(gca,[0 30]);set(gca,'XTick',bounds);
            print('-dpdf', fullfile(plots_dir,sprintf('Velocity_incorrect_%s_%s_%s_%s',ageism, name_reclayer, block{ddd}, genotypes{1})),'-painters');
            close
            figure
            graph_TG_incorrect = shadedErrorBar([], nanmean(group_speed_fails{analyze_groups(2)}),frequ_analysis.sem(group_speed_fails{analyze_groups(2)}), 'r-',1);hold off
            ylim(gca,[0 30]);set(gca,'XTick',bounds);
            print('-dpdf', fullfile(plots_dir,sprintf('Velocity_incorrect_%s_%s_%s_%s',ageism, name_reclayer, block{ddd}, genotypes{2})),'-painters');
            
            %%
            close
            figure
            graph_WT_correct = shadedErrorBar([], nanmean(group_speed{analyze_groups(1)},1),frequ_analysis.sem(group_speed{analyze_groups(1)}), 'k-',1);hold on
            graph_TG_correct = shadedErrorBar([], nanmean(group_speed{analyze_groups(2)},1),frequ_analysis.sem(group_speed{analyze_groups(2)}), 'g-',1);hold off            
            ylim(gca,[0 30]);
            set(gca,'XTick',bounds); %,'XTickLabel',{'return','delay','stem','choice','reward'});
            print('-dpdf', fullfile(plots_dir,sprintf('Velocity_%s_%s_%s',ageism, name_reclayer, block{ddd})),'-painters');
            close
            figure
            graph_WT_correct = plot(nanmean(group_speed{analyze_groups(1)},1), 'k-');hold on
            graph_TG_correct = plot(nanmean(group_speed{analyze_groups(2)},1), 'Color', [0 0.6 0]);hold off            
            ylim(gca,[0 30]);
            set(gca,'XTick',bounds); %,'XTickLabel',{'return','delay','stem','choice','reward'});
            print('-dpdf', fullfile(plots_dir,sprintf('Velocity_%s_%s_%s_line',ageism, name_reclayer, block{ddd})),'-painters');


            %% To plot power/frequ of beh and sleep1 of WT and TG
            line_color=  {[0 0 0];[0 0.52 0.25];[0 0 0];[0 0.52 0.25];[0 0 0];[0 0.2 0.9];};
            figure
            fig_WT_beh = shadedErrorBar(f,mean(power_frequ_beh{analyze_groups(1)},3)',frequ_analysis.sem(squeeze(power_frequ_beh{analyze_groups(1)})'),{'--','LineWidth',1 ,'Color',line_color{analyze_groups(1)}},1); hold on
            fig_WT_sleep = shadedErrorBar(f,mean(power_frequ_sleep{analyze_groups(1)},3)',frequ_analysis.sem(squeeze(power_frequ_sleep{analyze_groups(1)})'),{'LineWidth',1 ,'Color',line_color{analyze_groups(1)}},1); hold on
            fig_TG_beh = shadedErrorBar(f,mean(power_frequ_beh{analyze_groups(2)},3)',frequ_analysis.sem(squeeze(power_frequ_beh{analyze_groups(2)})'),{'--','LineWidth',1 ,'Color',line_color{analyze_groups(2)}},1);  hold on
            fig_TG_sleep = shadedErrorBar(f,mean(power_frequ_sleep{analyze_groups(2)},3)',frequ_analysis.sem(squeeze(power_frequ_sleep{analyze_groups(2)})'),{'LineWidth',1 ,'Color',line_color{analyze_groups(2)}},1);  hold off

            title(sprintf('power/frequ of beh and sleep of WT (n=%s) and TG (n=%s) %s', num2str(size(power_frequ_beh{analyze_groups(1)},3)),num2str(size(power_frequ_beh{analyze_groups(2)},3)), block{ddd}));hold off
            set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
            xlabel('Freq (Hz)'); ylabel('Raw Power');%ylim(gca,[0.2 2.5]);
            print('-dpdf', fullfile(plots_dir,sprintf('PowerFreq_WT_TG_%s', block{ddd})), '-painters')
            print('-dpng', fullfile(plots_dir,sprintf('PowerFreq_WT_TG_%s', block{ddd})), '-painters')
            
            %% Searching for error in power/frequ analysis
            line_color=  {[0 0 0];[0.0 0.9 0.5];[0 0 0];[0.0 0.9 0.5];[0 0 0];[0.2 0.2 0.74]};
            legendlabel= {'return','delay','stem','choice','reward','centr_arm'};
            for area = 1:6
                figure
                for gg = analyze_groups
                    for i = 1:size(power_freq_group_area{area,gg},3)
                        pwfq =  power_freq_group_area{area,gg};
                        plot(f, pwfq(:,:,i),'LineWidth',1 ,'Color',(line_color{gg}));hold on
                        %plot(f, mean(power_freq_group_area_incorr{area,gg},3),'--','LineWidth',1 ,'Color',(line_color{area}));
                    end
                    clear pwfq
                end
                title(sprintf('All Groups -- SleepNormalized -- %s -- %s', block{ddd},legendlabel{area}));hold off
                set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                xlabel('Freq (Hz)'); ylabel('Power normalized to sleep');ylim(gca,[0.2 2]);
                print('-dpdf', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_AllAnimals_%s_log_parsed_%s', block{ddd},legendlabel{area})))
                print('-dpng', fullfile(plots_dir,sprintf('PowerFreq_corrincorr_AllAnimals_%s_log_parsed_%s', block{ddd},legendlabel{area})))
                
            end
            
        end
    end
    
end


rmpath('Tools');