function ptempl = parsingtemplate(pmode)
%PARSINGTEMPLATE Define zones used to assign zone index to path data
%
% ptempl = PARSINGTEMPLATE(pmode)
%    pmode can be 'fig8:rdscr' (return-delay-stem-choice-reward)
%                 'fig8mouse:rdscr' (return-delay-stem-choice-reward)

% split = strsplit(pmode, ':');
% mazetype = split{1};
% pmode = split{2};

mazetype = 'fig8mouse';
pmode = 'rdscr';

switch mazetype
	case {'fig8rat'}
		if strcmpi(pmode, 'rdscr') % return-delay-stem-choice-reward
			ptempl(1).x = [-70 -70 -25 -25 -70 NaN 25 25 70 70 25];    % return left and right
			ptempl(1).y = [-90 60 60 -90 -90 NaN -90 60 60 -90 -90]; % return left and right
			ptempl(1).zone = 'return';
			
			ptempl(2).x = [-10 -10 10 10 -10]; % delay
			ptempl(2).y = -[90 45 45 90 90];   % delay
			ptempl(2).zone = 'delay';
			
			ptempl(3).x = [-20 -20 20 20 -20]; % stem
			ptempl(3).y = [-45 60 60 -45 -45]; % stem
			ptempl(3).zone = 'stem';
			
			ptempl(4).x = [-10 -10 10 10 -10]; % choice
			ptempl(4).y = [30 90 90 30 30];    % choice
			ptempl(4).zone = 'choice';
			
			ptempl(5).x = [-60 -60 -10 -10 -60 NaN 10 10 60 60 10]; % reward
			ptempl(5).y = [60 90 90 60 60 NaN 60 90 90 60 60];      % reward
			ptempl(5).zone = 'reward';
		end
		
		for i = 1:length(ptempl)
			[ptempl(i).x, ptempl(i).y] = poly2cw(ptempl(i).x, -ptempl(i).y);
		end
		
end
%% Fix mouse maze 22/10/2016 MH & SVS = <3
switch mazetype
	case 'fig8mouse'
        		if strcmpi(pmode, 'rdscr') % return-delay-stem-choice-reward
			ptempl(1).x = [-35 -35 -12 -12 -35 NaN 12 12 35 35 12 NaN -12 -12 12 12 -12];    % return left and right
			ptempl(1).y = [-50 25 25 -50 -50 NaN -50 25 25 -50 -50 NaN -50 -33 -33 -50 -50]; % return left and right
			ptempl(1).zone = 'return';
			
			ptempl(2).x = [-10 -10 10 10 -10]; % delay
			ptempl(2).y = -[33 15 15 33 33];   % delay
			ptempl(2).zone = 'delay';
			
			ptempl(3).x = [-10 -10 10 10 -10]; % stem
			ptempl(3).y = [-15 20 20 -15 -15]; % stem
			ptempl(3).zone = 'stem';
			
			ptempl(4).x = [-10 -10 10 10 -10]; % choice
			ptempl(4).y = [20 50 50 20 20];    % choice
			ptempl(4).zone = 'choice';
			
			ptempl(5).x = [-35 -35 -10 -10 -30 NaN 10 10 35 35 10]; % reward
			ptempl(5).y = [25 50 50 25 25 NaN 25 50 50 25 25];      % reward
			ptempl(5).zone = 'reward';
                end
        for i = 1:length(ptempl)
			[ptempl(i).x, ptempl(i).y] = poly2cw(ptempl(i).x, -ptempl(i).y);
        end
        %{
        % If you just want to resize
		for i = 1:length(ptempl)
			[ptempl(i).x, ptempl(i).y] = deal(ptempl(i).x*.5, ptempl(i).y*.5);
		end
        %}
        
        %% To look at this use:
            %{
                figure,
            for i = 1:5
            plot(ptempl(i).x,ptempl(i).y), hold on
            end
            plot(pathData{1,1}(3).x,pathData{1,1}(3).y)
            hold off
            %}
    
end