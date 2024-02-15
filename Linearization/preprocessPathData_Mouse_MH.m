% Modified MH 9/22/2016
addpath(genpath(fullfile(pwd, 'functions')));


%% Define here which animal / days to run
process_data = [1450] 

load('.\Figure8DataOrganization\sessionInfo.mat');
All_sessInfo = sessInfo;

for i = process_data
    
clearvars -except i process_data All_sessInfo
sessInfo = All_sessInfo(i); 
        
bInd = 1; % Starting index number; dummy variable for debuggin -- should be set to 1 otherwise
if bInd ~= 1
	warning('In DEBUG mode...')
end


overwrite = true;

parsemode = 'fig8mouse:rdscr';

regToUse = {'center','choice','base', 'return', 'reward', 'choiceEnd'}; % used to be just first 3: changed on @date 9/18/15 by @author: Sia
trialSuccess = 'correct';
trialDirection = 'any';

switch trialSuccess
    case 'any'
        successInclude = {'success',[true false]};
    case 'correct'
        successInclude = {'success',true};
    case 'incorrect'
        successInclude = {'success',false};
end
switch trialDirection
    case 'any'
        dirInclude = {'direction',{'L','R'}};
    case 'left'
        dirInclude = {'direction',{'L'}};
    case 'right'
        dirInclude = {'direction',{'R'}};
end

regLocs.center = {'A25'};
regLocs.choice = {'N5'};    %
regLocs.base = {'A56','A45'};
regLocs.return = {'A16','N1','A12','A23','N3','A34'};
regLocs.reward = {'N6','N4'};
regLocs.choiceEnd = {'N2'}; %

allRegions = fields(regLocs);
locsToUse = [];
for r = 1:length(regToUse)
   locsToUse = cat(2,locsToUse,regLocs.(regToUse{r})); 
end

epochDirs = {sessInfo.sessDirs}';  
%epochDirs{1,1}{1,1} = 'd0';  %%added FOR Debugging only
mainDirs = {sessInfo.mainDir}';
TFile = {sessInfo.tList}';

locInfo = cell(length(mainDirs), 1);
pathData = cell(length(mainDirs), 1);
pathDataIdeal = cell(length(mainDirs), 1);
pathDataLin = cell(length(mainDirs), 1);
trialInfo = cell(length(mainDirs), 1);
parseInfo = cell(length(mainDirs), 1);
didntWork = [];
for dirNum = bInd:length(mainDirs)
		
        sessDirs = cellfun(@(x) fullfile(mainDirs{dirNum}, x), epochDirs{dirNum}(:), 'UniformOutput', false);
        	disp(sessDirs);
        disp('-> loading path:'); disp(mainDirs{dirNum});
          pathData{dirNum} = loadpaths_indata(mainDirs{dirNum});
        
       
		disp('-> Path loaded');
        %%
        	disp('-> Idealizing Path');
       [pathDataIdeal{dirNum}, locInfo{dirNum}, maze] = idealizepath(pathData{dirNum}, 'fig8mouse');
        	disp('-> Done');
		trialInfo{dirNum} = trialsFromLoc(locInfo{dirNum});
        
		pathDataLin{dirNum} = linearizepath(pathDataIdeal{dirNum},trialInfo{dirNum});
		parseInfo{dirNum} = parsepath(pathDataIdeal{dirNum},parsingtemplate(parsemode));
		parseInfo{dirNum} = parseinfo2trial(parseInfo{dirNum}, trialInfo{dirNum});
        
        figure
        plot(pathDataIdeal{1,1}(2).x,pathDataIdeal{1,1}(2).y)
        
end
%%
% separate by subtrial (block) and save in respective folder:
% locInfo, pathData, pathDataIdeal, pathDataLin, trialInfo
rememberWhoDidntGetSaved = [];
for dirNum = setdiff(bInd:length(mainDirs), didntWork)
	sessDirs = cellfun(@(x) fullfile(mainDirs{dirNum}, x), epochDirs{dirNum}, 'UniformOutput', false);

	for s = 1:length(sessDirs)
		li = locInfo{dirNum}(s);
		pd = pathData{dirNum}(s);
		pdi = pathDataIdeal{dirNum}(s);
		pdl = pathDataLin{dirNum}(s);
		ti = trialInfo{dirNum}(s);
		
		pri.inds = parseInfo{dirNum}{s};
		pri.tInt = ti.tInt;
		pri.direction = ti.direction;
		pri.success = ti.success;
		
		if ~exist(fullfile(sessDirs{s}, 'locInfo.mat'), 'file') || overwrite
			save(fullfile(sessDirs{s}, 'locInfo.mat'),'-struct','li');
		else
			rememberWhoDidntGetSaved{end+1} = fullfile(sessDirs{s}, 'locInfo.mat');
		end
		if ~exist(fullfile(sessDirs{s}, 'pathData.mat'), 'file') || overwrite
			save(fullfile(sessDirs{s}, 'pathData.mat'),'-struct','pd');
		else
			rememberWhoDidntGetSaved{end+1} = fullfile(sessDirs{s}, 'pathData.mat');
		end
		if ~exist(fullfile(sessDirs{s}, 'pathDataIdeal.mat'), 'file') || overwrite
			save(fullfile(sessDirs{s}, 'pathDataIdeal.mat'),'-struct','pdi');
		else
			rememberWhoDidntGetSaved{end+1} = fullfile(sessDirs{s}, 'pathDataIdeal.mat');
		end
		if ~exist(fullfile(sessDirs{s}, 'pathDataLinear.mat'), 'file') || overwrite
			save(fullfile(sessDirs{s}, 'pathDataLinear.mat'),'-struct','pdl');
		else
			rememberWhoDidntGetSaved{end+1} = fullfile(sessDirs{s}, 'pathDataLinear.mat');
		end
		if ~exist(fullfile(sessDirs{s}, 'trialInfo.mat'), 'file') || overwrite
			save(fullfile(sessDirs{s}, 'trialInfo.mat'),'-struct','ti');
		else
			rememberWhoDidntGetSaved{end+1} = fullfile(sessDirs{s}, 'trialInfo.mat');
		end
		if ~exist(fullfile(sessDirs{s}, 'parsingInfo.mat'), 'file') || overwrite
			save(fullfile(sessDirs{s}, 'parsingInfo.mat'),'-struct','pri');
		else
			rememberWhoDidntGetSaved{end+1} = fullfile(sessDirs{s}, 'trialInfo.mat');
		end
	end
end

  %  catch %try in line 19
  %      errordlg(sprintf('Could not process i %d Animal %d Day %d', i,sessInfo.animal,sessInfo.day));
  %  end %try in line 19
    
end         % end of cycling through days

rmpath(genpath(fullfile(pwd, 'functions')));