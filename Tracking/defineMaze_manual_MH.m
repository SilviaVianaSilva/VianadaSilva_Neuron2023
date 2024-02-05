% Define Maze
% Standalone script to define the maze manually when the path is rotated
% and maze arm-lengths mismatch. Circulate through new data and save file
% for a recording day in processedData.
%   Matthias Haberl 09/15/2016

clearvars -except previous_i
%addpath(genpath('C:\Users\Silvia\Desktop\Newest MClust (microvolts, history logging)'));
addpath(genpath('C:\Users\Silvia\Dropbox\SilviaProjectCode\MousePath'));
load('C:\Users\Silvia\Dropbox\SilviaProjectCode\+Figure8DataOrganization\figure8Data.mat');

for i = [1450]  
    % 410 453 304 314 614 873     1623 1629 1631 1635 2259 2250 2554 
    
  
sessDirs = path{i}; savefolder = fullfile(sessDirs,'processedData');mkdir(savefolder); 

sessNames = bSess{i};       
    for ii = 1%:n_sessions

    disp(strcat('Start reading the input data for: ',fullfile(sessDirs, sessNames{ii})));

    %% Get position data
    videoFile = fullfile(sessDirs,sessNames{ii},'VT1.Nvt');
   
    [x,y,t,angle,n_samples(ii),n_trackingPoints(ii),n_NaNPoints(ii)] ...
        = writeVelocity.getPreprocessedPositionData(videoFile,diode(i));              

    
    %% How do you define your first values
    if exist(fullfile(path{i},'processedData', 'maze_coordinates.mat'), 'file')==2
        load(fullfile(path{i},'processedData', 'maze_coordinates.mat'));
    elseif exist('previous_i', 'var')
        try
           filename = fullfile(path{previous_i},'processedData', 'maze_coordinates.mat');   
           load(filename);
        end 
    end    
    if exist('coord', 'var')
        SW = coord.SW; SE= coord.SE; NE= coord.NE; NW = coord.NW;
    else    
            % Find border values for path and box
            maxX = max(x);
            minX = min(x);
            maxY = max(y);
            minY = min(y);

            % Set the corners of the reference box
            NE = [maxX, maxY]
            NW = [minX, maxY]
            SW = [minX, minY]
            SE = [maxX, minY]    
            clear coord.line
            linecoordinates = [mean([NE;NW]);mean([SE;SW])];
            coord.line = [linecoordinates(:,1),linecoordinates(:,2)]
    end
f = figure;
 plot(x,y);
 axis equal
 title(num2str(i))
 max(x);
    h = impoly(gca, [SW; SE; NE; NW])
    setColor(h,'red');    
    %addNewPositionCallback(h,@(p) title(mat2str(p,3)));
    fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn);

    ll = imline(gca,coord.line); % 
    setColor(ll,[1 0 1]);
    %setString(ll,'center');
 %{
    pp = impoint(gca,100,200)
    % Update position in title using newPositionCallback
    addNewPositionCallback(pp,@(pp) title(sprintf('(%1.0f,%1.0f)',pp(1),pp(2))));
    % Construct boundary constraint function
    fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
    % Enforce boundary constraint function using setPositionConstraintFcn
    setPositionConstraintFcn(pp,fcn);
    setColor(pp,'r');
    setString(pp,'center');
    %}
    
continue_button = uicontrol('Position', [1 1 200 40], 'String', 'Continue', ...
                      'Callback', 'uiresume(gcbf)');
predict_button = uicontrol('Position', [250 1 200 40], 'String', 'Preview', ...
                      'Callback', 'writeVelocity.Maze_Preview_MH');
%closeall_button = uicontrol('Position', [500 1 100 40], 'String', 'Close All', ...
%                      'Callback', 'close all');            
disp('Waiting for user to select maze coordinates');
uiwait(f);
pos = getPosition(h)
coord.SW = [pos(1,1), pos(1,2)];
coord.SE = [pos(2,1), pos(2,2)];
coord.NE = [pos(3,1), pos(3,2)];
coord.NW = [pos(4,1), pos(4,2)];
%disp(sprintf('Coordinates %d', h)); 
coord.line = getPosition(ll);

coord.center_rot = atand((max(coord.line(:,1))-min(coord.line(:,1)))/(max(coord.line(:,2)-min(coord.line(:,2)))))

coord.west_arm_rot = atand((coord.NW(1)-coord.SW(1))/(coord.NW(2)-coord.SW(2)))
coord.east_arm_rot = atand((coord.NE(1)-coord.SE(1))/(coord.NE(2)-coord.SE(2)))



coord.center_rot = mean([coord.center_rot, coord.west_arm_rot, coord.east_arm_rot]);
coord.theta_deg = (coord.center_rot);
coord.theta_rad = degtorad(coord.theta_deg); 
close all %close(f);
filename = fullfile(savefolder, 'maze_coordinates.mat');
save(filename, 'coord') 

    end
previous_i = i;    
end

%rmpath(genpath('C:\Users\Silvia\Desktop\Newest MClust (microvolts, history logging)'));
rmpath(genpath('C:\Users\Silvia\Dropbox\SilviaProjectCode\MousePath'));

