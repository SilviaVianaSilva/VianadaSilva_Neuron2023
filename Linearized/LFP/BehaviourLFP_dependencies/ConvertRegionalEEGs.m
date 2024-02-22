% Convert Regional EEGs and Normalize to quiet S1 of the same day of this animal
% therefore needs to load 'quiet_SPG' from mainDir/s1/LFP/Spectogram_channel%d.mat
% => needs to 
%
%
% run ExtractRegionalEEG_forSpectAnalysis before
% run LFP_Analysis_Sleep before
%
%  Input:
% - EEG metafile(to know if SPG has been extracted) and animalnumbers to analyse
% - behaviour SPG
% - quiet_SPG
%
% run Analyze_AverageSpectograms_sleepNorm afterwards
%
%---------------------------------------------------------
% Matthias Haberl, 2017
%---------------------------------------------------------

%% Load Metafile of all Regional EEGs extracted so far
load('metafile_RegionalEEGs.mat');


plot_per_day = 0; %%Turn on plot per day here 
%n_total = numel(animals_analysed);

load('+Figure8DataOrganization\sessionInfo.mat');
All_sessInfo = sessInfo; clear sessInfo

%% Load updated binning per time
frequ_analysis.load_TimeBins

%% Initialize
failed_iii = [];
timer_prog = clock;
block = {'d0','d2','d10'};
all_rois = {'return';'delay';'stem';'choice';'reward'};

d0_info = struct2table(metainfo.d0);

for ddd = 1:3 % delay d0 d2 d10, never to be directly compared, also called block{}
    fprintf('---> Processing %s <---\n', block{ddd});
    binsizes =  opt.binninByTime.(block{ddd});
    for aaa = animals_analysed  % Animals to process
        animal_spect = [];  animal_speed = [];
        animal_speed_p_region = [];
        animal_sleepNormalized = []; animal_fails_sleepNormalized = [];
       
        animal_fail_spect = [];  animal_speed_fail = [];
        animal_quiet_power_spect = []; 
        animal_active_power_spect=[];
         
        timer_anim = tic;
        i_list = find((animal_no==aaa) & (group_A | group_B))  % groups generated in datadef2groups, like this recording layers aren't mixed etc.
        is_age = All_sessInfo(i_list(1)).age; %Check the age of this animal
        reclayer = All_sessInfo(i_list(1)).tListLoc;
        for iii = i_list
            sessInfo = All_sessInfo(iii);
            
            
            if ~isempty(metainfo.(block{ddd})(iii).filename) %doublecheck that this file was processed
                
                %try
                fprintf('\n -- Starting to process i number: %d --\n', iii);
                close all
                clearvars spect day_spect adjusted_spects_corr adjusted_spects_incorr unsuccessf_day_spect speed speed_fail day_speed fail_day_speed
                clearvars  sleepNormalized_day_spect sleepNormalized_fails_day_spect
                %% Load spectogram
                if exist(metainfo.(block{ddd})(iii).filename, 'file')==2
                load(metainfo.(block{ddd})(iii).filename);                
                else
                    errordlg(sprintf('Skipping missing file: %s' ,metainfo.(block{ddd})(iii).filename));
                    continue
                end
                
                %% Loop through trials and regions
                suc_num=0;fail_num=0;
                f=frequ{1,1};
                adjusted_spects_corr=zeros(size(f,2), sum(binsizes),sum(successful));
                fail_lapse = (numel(successful)-sum(successful));
                
                if fail_lapse > 0; adjusted_spects_incorr=nan(size(f,2), sum(binsizes),fail_lapse); end
                for trial = 1:size(spect,1)
                    if successful(trial)  %use only succesfull trials
                        suc_num=suc_num+1;
                        adjusted_spects_corr(:,:,suc_num) = [imresize(spect{trial,1},[size(f,2), binsizes(1)]),imresize(spect{trial,2},[size(f,2), binsizes(2)]),imresize(spect{trial,3},[size(f,2), binsizes(3)]),imresize(spect{trial,4},[size(f,2), binsizes(4)]),imresize(spect{trial,5},[size(f,2), binsizes(5)])];
                        speed(:,:,suc_num) = [imresize((velocity{trial,1})',[1,binsizes(1)]),imresize((velocity{trial,2})',[1,binsizes(2)]),imresize((velocity{trial,3})',[1,binsizes(3)]),imresize((velocity{trial,4})',[1,binsizes(4)]),imresize((velocity{trial,5})',[1,binsizes(5)])];
                        speed_p_region(:,:,suc_num) = [nanmean(velocity{trial,1}),nanmean(velocity{trial,2}),nanmean(velocity{trial,3}),nanmean(velocity{trial,4}),nanmean(velocity{trial,5})];
                    elseif successful(trial)==0  %use only unsuccesfull trials  
                        fail_num=fail_num+1;
                        adjusted_spects_incorr(:,:,fail_num) = [imresize(spect{trial,1},[size(f,2), binsizes(1)]),imresize(spect{trial,2},[size(f,2), binsizes(2)]),imresize(spect{trial,3},[size(f,2), binsizes(3)]),imresize(spect{trial,4},[size(f,2), binsizes(4)]),imresize(spect{trial,5},[size(f,2), binsizes(5)])];                   
                        speed_fail(:,:,fail_num) = [imresize((velocity{trial,1})',[1,binsizes(1)]),imresize((velocity{trial,2})',[1,binsizes(2)]),imresize((velocity{trial,3})',[1,binsizes(3)]),imresize((velocity{trial,4})',[1,binsizes(4)]),imresize((velocity{trial,5})',[1,binsizes(5)])];                   
                    end

                end   
                if fail_num == 0
                    fprintf('No wrong trials\n'); 
                end
                day_spect = mean(adjusted_spects_corr,3);
                day_speed = mean(speed,3);
                day_speed_p_region = mean(speed_p_region,3);
                 if fail_num > 0; unsuccessf_day_spect =  nanmean(adjusted_spects_incorr,3); fail_day_speed = nanmean(speed_fail,3); 
                 else
                     unsuccessf_day_spect = NaN(size(day_spect));fail_day_speed = NaN(size(day_speed));
                 end

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
            analysis_dir = ['.\LFP_Analysis\Sleep_MatFiles_PowerOverFreq\SleepPowerFreq_', Ageism,expressionlayer]; if ~exist(analysis_dir,'dir'), mkdir(analysis_dir); end
            powfrq_analysisfile = fullfile(analysis_dir,sprintf('animal_%d_i%d_s1.mat',sessInfo.animal,iii))
            
            load(powfrq_analysisfile, 'power_spect','quiet_power_spect','active_power_spect','quiet_power_std', 'active_power_std','f','sessInfo')
            
            %% Regular z-scored spect
            zsored_day_spect = (day_spect - repmat(mean(day_spect,2),1,size(day_spect,2))) ./ (repmat(std(day_spect,0,2),1,size(day_spect,2)));
            
            %% Dividing by S1
            sleepNormalized_day_spect = (day_spect ./ repmat(quiet_power_spect,1,size(day_spect,2)));
            
            %% Plot per day
            if plot_per_day    
                %outdir = fullfile('C:\Users\Silvia\Desktop\Plots\Spectogram_Methods' , sprintf('%s_day_%s', block{ddd}, num2str(iii)));mkdir(outdir)
                bounds = frequ_analysis.bins2bounds(binsizes);
                outdir = fullfile(sessInfo.mainDir,block{ddd},'Spectograms');mkdir(outdir);                                
                figure
                h = pcolor([0:size(day_spect,2)-1],f,  day_spect);
                h.EdgeColor = 'none';    
                set(gca,'YScale','log','XTick',bounds, 'YTick', [2, 4, 6, 8,16, 30, 100,200,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');axis xy, colorbar, title('single day spect')           
                filename = fullfile(outdir , sprintf('day_spect_%s_%s', num2str(iii),block{ddd}));
                print('-dpng', filename,'-opengl') %-painters :painters is vectorgraphics but slow and too big
                saveas(gcf,filename,'epsc')
                close all

                %% Regular z-scored spect (for comparison)
                figure
                h = pcolor([0:size(zsored_day_spect,2)-1],f,  zsored_day_spect);
                h.EdgeColor = 'none';        
                set(gca,'YScale','log','XTick',bounds, 'YTick', [2, 4, 6, 8,16, 30, 100,200,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                axis xy,colorbar, title('z-scored - single day spect')  
                filename = fullfile(outdir , sprintf('z-scored_%s_%s', num2str(iii),block{ddd}));
                print('-dpng', filename,'-opengl') %-painters :painters is vectorgraphics but slow and too big
                saveas(gcf,filename,'epsc')
                close all

                %% Dividing by S1 - Normalized to S1
                figure                 
                h = pcolor([0:size(sleepNormalized_day_spect,2)-1],f,  sleepNormalized_day_spect);
                h.EdgeColor = 'none';        
                set(gca,'YScale','log','XTick',bounds, 'YTick', [2, 4, 6, 8,16, 30, 100,200,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
                axis xy, colorbar, title('dividing by mean of sleep per frequ - single day spect')
                filename = fullfile(outdir , sprintf('normalizedtoS1_division_%s_%s', num2str(iii),block{ddd}));
                print('-dpng', filename,'-opengl') 
                saveas(gcf,filename,'epsc')
                close all  

            end % end plot_per_day
                
                % Unsuccesfull trials, SleepNormalized LFP
                if fail_num > 0; 
                    sleepNormalized_fails_day_spect = (unsuccessf_day_spect ./ repmat(quiet_power_spect,1,size(day_spect,2)));

                else
                     sleepNormalized_fails_day_spect = NaN(size(day_spect));
                end
                
            else
                errordlg(sprintf('Behav EEG of i: %s was not processed', num2str(iii)));
                continue
            end
            animal_spect  = cat(3,animal_spect,day_spect);
            animal_fail_spect = cat(3,animal_fail_spect,unsuccessf_day_spect);
            
            animal_speed = cat(3,animal_speed,day_speed);
            animal_speed_fail = cat(3,animal_speed_fail,fail_day_speed);
            
            animal_speed_p_region = cat(3,animal_speed_p_region,day_speed_p_region);
            
            animal_sleepNormalized  = cat(3, animal_sleepNormalized, sleepNormalized_day_spect);
            animal_fails_sleepNormalized  = cat(3, animal_fails_sleepNormalized, sleepNormalized_fails_day_spect);
            
            animal_quiet_power_spect = cat(3,animal_quiet_power_spect, quiet_power_spect);
            animal_active_power_spect = cat(3,animal_active_power_spect,active_power_spect);
            
            clear quiet_power_spect active_power_spect quiet_power_std active_power_std
        end
        close all

        
        %Control for each animal 
        outfolder = ('.\Plots\Aver_spect_sleepNormBehav\power_frequ'); mkdir(outfolder);
        filename = fullfile(outfolder,sprintf('Animal_%d_power_frequ_recLay_%s_%s',aaa,reclayer,block{ddd}));
        frequ_analysis.plot_Behav_quietSl_MovSl(filename,f,mean(animal_spect,3),mean(animal_sleepNormalized,3),mean(animal_quiet_power_spect,3),mean(animal_active_power_spect,3));        
            
        outfolder = ('.\LFP_Analysis\AnimalSpects');mkdir(outfolder);
        filename = fullfile(outfolder,sprintf('Animal_%d_adjusted_spect_behav_recLay_%s_%s',aaa,reclayer,block{ddd}));
        save(filename,'animal_spect','f','animal_sleepNormalized','animal_fail_spect','animal_fails_sleepNormalized','animal_speed_p_region','animal_speed', 'animal_speed_fail','animal_quiet_power_spect');
        
        figure
        plot(nanmean(animal_speed,3),'k','LineWidth', 1), hold on
        ylabel('Speed');xlabel('Location');
        title('Speed');
        filename = fullfile(outfolder,sprintf('Animal_%d_AverSpeed_scaled2Spectogr_recLay_%s_%s',aaa,reclayer,block{ddd}));
        print('-dpdf', filename)
        
        
       
        
    end
    
end

elapsed_min = round(etime(clock,timer_prog)/60);
msgbox(sprintf('Total runtime: %d min\n',elapsed_min),'Runtime');