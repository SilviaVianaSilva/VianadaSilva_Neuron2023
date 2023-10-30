%% Plot Power over frequency

clear
addpath('Tools'); 

analyze = 2;  % 1: Analyze old animals, 2: Analyze Young Animals




load('+Figure8DataOrganization\sessionInfo.mat');
plots_BaseDir = ('.\Plots\sleepLFP'); if ~exist(plots_BaseDir,'dir'), mkdir(plots_BaseDir),end
All_sessInfo = sessInfo; clear sessInfo
animal_numbers = unique([All_sessInfo.animal]);
n_total = numel(animal_numbers);

All_animal_data =struct2dataset(All_sessInfo);
WT_assembled_power_spect_s1=[];WT_assembled_power_spect_s2=[];
TG_assembled_power_spect_s1=[];TG_assembled_power_spect_s2=[];
%% Go through animals
for aaa = animal_numbers
    clearvars power_spect is_I i_list
    assembled_power_spect_s1 = []; assembled_power_spect_s2 = [];
    for t = 1:size(All_sessInfo,1)
     is_I(t)=ismember(aaa,All_sessInfo(t).animal); %pick days of this animal
    end
    i_list=find(is_I);
    is_age = All_sessInfo(i_list(1)).age; %Check the age of this animal   
    
    if analyze == 1
        if is_age<10, fprintf('\nSkipping Young Animal %d, Age:%d ',aaa ,is_age) , continue, end  %skip if wrong category
        Ageism = 'old';
        plots_Dir = fullfile(plots_BaseDir,Ageism);
    elseif analyze == 2
        if is_age>10,fprintf('\nSkipping Old Animal %d, Age:%d ',aaa ,is_age), continue, end
        %elseif sessInfo.age<10
        Ageism = 'young';
        plots_Dir = fullfile(plots_BaseDir,Ageism);
    end
    if ~exist(plots_Dir,'dir'), mkdir(plots_Dir),end
    analysis_dir = ['.\MATLAB\SleepPowerFreq_', Ageism];
    
    animals_genotype = All_sessInfo(i_list(1)).genotype;
    
    %% Go through i-list of this animal
    for iii = i_list
        filename = fullfile(analysis_dir,sprintf('animal_%d_i%d_%s.mat',aaa,iii,'s1'));
        if exist(filename,'file')
            load(filename, 'power_spect','f','sessInfo','block')
            assembled_power_spect_s1 = [assembled_power_spect_s1;power_spect];
        else
            fprintf('\nCould not find: %s', filename)
        end
        clear power_spect      
        filename = fullfile(analysis_dir,sprintf('animal_%d_i%d_%s.mat',aaa,iii,'s2'))
        if exist(filename,'file')
        load(filename, 'power_spect','f','sessInfo','block')
        assembled_power_spect_s2 = [assembled_power_spect_s2;power_spect];
        else
            fprintf('\nCould not find: %s', filename)
        end        
    end
    
    if isempty(assembled_power_spect_s1),fprintf('\n No Spectra found from Animal %d ',aaa), continue, end
    %% Plot curve for this animal
    figure
    subplot(2,1,1)
    %scatter(bounds,[0 0 0 0 0 0]), hold on
    s1_plot = shadedErrorBar(f,mean(assembled_power_spect_s1,1),sem(assembled_power_spect_s1),{'k','Linewidth',2},1),hold on
    s2_plot = shadedErrorBar(f,mean(assembled_power_spect_s2,1),sem(assembled_power_spect_s2),{'g','Linewidth',2},1),hold off
    title(sprintf('S1 and S2 animal no. - %d',aaa));
    xlim(gca,[4 300]);
    ylim(gca,[0 0.0012]);
    set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
    ylabel('Mean Power'); xlabel('Freq (Hz)');
    legend([s1_plot.mainLine,s2_plot.mainLine], {'S1','S2'},'Location','NorthEast');
    %% Fitting log curve
    %data = (mean(assembled_power_spect_s1,1));
    %f2 = fit(f',data','exp1')
    %plot(f2,f',data')
    %set(gca, 'yscale', 'log', 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
   
    %axis tight
    filename = fullfile(plots_Dir,sprintf('Power_frequ_s1s2_Animal%02d', aaa));
    savefig(filename)
    print('-dpdf', filename)
    %close all
    
    %% Assemble this animal into all WT and TG animals
    if animals_genotype == 1
        WT_assembled_power_spect_s1 = [WT_assembled_power_spect_s1; mean(assembled_power_spect_s1,1)]
        WT_assembled_power_spect_s2 = [WT_assembled_power_spect_s2; mean(assembled_power_spect_s2,1)]
    elseif animals_genotype == 2
        TG_assembled_power_spect_s1 = [TG_assembled_power_spect_s1; mean(assembled_power_spect_s1,1)]
        TG_assembled_power_spect_s2 = [TG_assembled_power_spect_s2; mean(assembled_power_spect_s2,1)]
    else
        warndlg('Undefined genotype');
    end
end
    %% Plot all WT and TG animals
    
    figure
    subplot(2,1,1)
    %scatter(bounds,[0 0 0 0 0 0]), hold on
    WT_plot1 = shadedErrorBar(f,mean(WT_assembled_power_spect_s1,1),sem(WT_assembled_power_spect_s1),{'k','Linewidth',2},1),hold on
    TG_plot1 = shadedErrorBar(f,mean(TG_assembled_power_spect_s1,1),sem(TG_assembled_power_spect_s1),{'r','Linewidth',2},1),hold off
    title(sprintf('S1 - WT and TG animals - %s',Ageism));
    xlim(gca,[4 300]); ylim(gca,[0 0.0012]);
    set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
    ylabel('Mean Power'); xlabel('Freq (Hz)');
    legend([WT_plot1.mainLine,TG_plot1.mainLine], {'WT','TG'},'Location','NorthEast');
    axis tight
    subplot(2,1,2)
    %scatter(bounds,[0 0 0 0 0 0]), hold on
    WT_plot2 = shadedErrorBar(f,mean(WT_assembled_power_spect_s2,1),sem(WT_assembled_power_spect_s2),{'k','Linewidth',2},1),hold on
    TG_plot2 = shadedErrorBar(f,mean(TG_assembled_power_spect_s2,1),sem(TG_assembled_power_spect_s2),{'r','Linewidth',2},1),hold off
    title(sprintf('S2 - WT and TG animals - %s',Ageism));
    xlim(gca,[4 300]); ylim(gca,[0 0.0012]);
    set(gca, 'xscale', 'log', 'XTick', [4, 6, 8, 30, 100,300], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'TickDir', 'out');
    ylabel('Mean Power'); xlabel('Freq (Hz)');
    legend([WT_plot2.mainLine,TG_plot2.mainLine], {'WT','TG'},'Location','NorthEast');
    axis tight    
    
    filename = fullfile(plots_BaseDir,sprintf('Power_frequ_s1s2_%s', Ageism));
    savefig(filename)
    print('-dpdf', filename)
    close all