function [indata, scale] = getPreprocessedIndata(sessDirs,sessNames,doublediode,target,onMaze,varargin)
extract_varargin;
indata = [];
n_sessions = length(sessNames);
n_trackingPoints = zeros(n_sessions,1);
n_NaNPoints = zeros(n_sessions,1);
n_samples = zeros(n_sessions,1);
for ii = 1:n_sessions
    %% Display progress
    fprintf('Start reading the input data for: \n',fullfile(sessDirs, sessNames{ii}));
   
    %% Get position data
    videoFile = fullfile(sessDirs,sessNames{ii},'VT1.Nvt');
   
    if ~exist(videoFile,'file')
    errordlg(sprintf('Missing file: %s',videoFile));
    end
    [x,y,t,angle,n_samples(ii),n_trackingPoints(ii),n_NaNPoints(ii)] ...
        = writeVelocity.getPreprocessedPositionData(videoFile,doublediode);

        if onMaze == 1
                    
            save_dir = fullfile(videoFile(1:end-7),'TrackingImages');
            mkdir(save_dir);
            [x,y, centerX, centerY,maze,scale] = writeVelocity.rotate_maze_MH(x, y, sessDirs, target, save_dir);  
            
            scale
            centerX = 0; centerY = 0;
            
        elseif onMaze == 2 % For open field == 2
            save_dir = fullfile(videoFile(1:end-7),'TrackingImages');
            mkdir(save_dir);
            centerX = (min(x)+max(x))/2;
            centerY = (min(y)+max(y))/2;
            x = (x - centerX);
            y = (y - centerY);
            
            y_scaling = target(2) / (max(y)-  min(y));
            x_scaling = target(1) / (max(x)-  min(x));
            scale = mean([y_scaling, x_scaling]);
             x = x * scale ;
            y = y * scale;
        else % for sleep == 0
            
        %% Adjust sleep box data
        scale = scaleB;
        
        clear center*
        save_dir = fullfile(videoFile(1:end-7),'TrackingImages');
        mkdir(save_dir);
        figure
        plot(x,y,'k');title('Sleep Before processing'); axis equal 
        saveas(gcf, fullfile(save_dir,'BeforeProcessing.png')); 
            centerX = (min(x)+max(x))/2;
            centerY = (min(y)+max(y))/2;     
        x = (x - centerX);
        y = (y - centerY);
        close
                figure
        plot(x,y,'k');title('Sleep Centered'); axis equal 
        saveas(gcf, fullfile(save_dir,'Rotated_NotScaled.png')); 
        close
        x = x * scaleB ;
        y = y * scaleB;
    
                figure
        plot(x,y,'k');title('Sleep Centered and Scaled'); axis equal       
        saveas(gcf, fullfile(save_dir,'Rotated_AndScaled.png')); 
        close

        maze.center = [mean(x) mean(y)];
            % x = (x - centerX) * scaleB(1);
            % y = (y - centerY) * scaleB(2);
            % scale = scaleB;
        end       

fh = figure();
    switch doublediode 
    case 1
    plot(x,y,'b');       
    case 2    
    plot(x,y,'r'); 
    case 3
    plot(x,y,'g');
    otherwise
    plot(x,y,'k');            
    end
    
axis equal    
sequ = 1; fname = fullfile(save_dir, sprintf('finalVelocity_Tracking.pdf'));  
print(fh,fname,'-dpdf');
close(fh);
%}



%-------------------------------------------------------
%% Put all data into the struct array
%-------------------------------------------------------
    
    indata = [indata; struct('x',x', 'y',y', 't',t', 'angle',angle)]; %,'maze',maze,'centerX',centerX,'centerY',centerY)];

    
end %sessions loop

%-------------------------------------------------------
%%addVelocity
%-------------------------------------------------------
[indata speed] = writeVelocity.addVelocityToIndata(indata);

%-------------------------------------------------------
%% Save indata
%-------------------------------------------------------
mkdir(strcat(sessDirs,'\processedData'));
cd([sessDirs,'\processedData']);
if onMaze == 1
    save('indataB','indata','n_samples','n_trackingPoints','n_NaNPoints','target');
elseif onMaze == 2
    save('indata_of1','indata','n_samples','n_trackingPoints','n_NaNPoints','target');
elseif onMaze == 0
    save('indataS','indata','n_samples','n_trackingPoints','n_NaNPoints','target');
end