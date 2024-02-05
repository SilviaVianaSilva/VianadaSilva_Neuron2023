function [x_rotated_scaled,y_rotated_scaled, x_center,y_center,maze,scaling_for_sleep] = rotate_maze_MH(x, y, sessDirs, target_maze, save_dir);  

% Rotates animal path around a central point. Central point and rotation
% angle are defined using defineMaze_manual_MH.
% The rotation angle and central point are calculated from the central axis.
% Then performs a scaling that is depending on the corner points defined in
% the user interface of defineMaze.

   
filename = fullfile(sessDirs,'processedData', 'maze_coordinates.mat');   
load(filename);

%% Rotate Maze
  
%-----------------------------------------------------
%% Get Angle; Set center
%-----------------------------------------------------

theta = coord.theta_rad;

% choose a point which will be the center of rotation
x_center = mean(coord.line(:,1));
y_center = mean(coord.line(:,2));

% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(x));

%angle = %rotation_matrix = repmat([
% define a 60 degree counter-clockwise rotation matrix %theta = pi/3;       % pi/3 radians = 60 degrees

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

%-----------------------------------------------------
%% Rotate running path
%-----------------------------------------------------
v = [x;y];          % create a matrix to apply
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
x_rotated = so(1,:);
y_rotated = so(2,:);

%% Taking this line out will let your maze be centered
% vo = so + center;   % shift again so the origin goes back to the desired center of rotation                     
% this can be done in one line as:
% vo = R*(v - center) + center

%-----------------------------------------------------
%% Now rotate poly
%-----------------------------------------------------

NE = [coord.NE(1); coord.NE(2)];NW = [coord.NW(1); coord.NW(2)];
SE = [coord.SE(1); coord.SE(2)];SW = [coord.SW(1); coord.SW(2)];
maze.NE_rot = R*(NE - [x_center; y_center]) %+ [x_center; y_center] %put either back to origin, or taking second part out, new center is 0
maze.NW_rot = R*(NW - [x_center; y_center]) %+ [x_center; y_center]
maze.SE_rot = R*(SE - [x_center; y_center]) %+ [x_center; y_center]
maze.SW_rot = R*(SW - [x_center; y_center]) %+ [x_center; y_center]
maze.xCenter = 0; maze.yCenter = 0;
% pick out the vectors of rotated x- and y-data

%-----------------------------------------------------
%% Make plot: Rotated not scaled
%-----------------------------------------------------
fh = figure()
%plot(x, y, 'b-', x_rotated, y_rotated, 'r-', x_center, y_center, 'bo');%
%hold on
scatter([maze.SW_rot(1), maze.SE_rot(1), maze.NW_rot(1), maze.NE_rot(1)], [maze.SW_rot(2), maze.SE_rot(2), maze.NW_rot(2), maze.NE_rot(2)], 'g')
hold on
%%use tocompare with and without rotation
plot(x_rotated, y_rotated, 'r-', 0, 0, 'bo');
axis equal
hold off
title('Rotated not scaled');
saveas(fh,fullfile(save_dir,'Rotated_NotScaled.png')); 
close

%-----------------------------------------------------
%% Scaling
%-----------------------------------------------------
target(1) = target_maze(1) / 2; %divide target_size by 2 to get size of one arm
target(2) = target_maze(2); %divide target_size by 2 to get size of one arm


%% 1: Scaling Arms proportionally in x and in y along the entire axis
positive_x = find(x_rotated>0);negative_x = find(x_rotated<0);
x_rotated_scaled(positive_x) = x_rotated(positive_x) .* (target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)])) ) ; x_scaling(1) =(target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)])) ) ;
x_rotated_scaled(negative_x) = -x_rotated(negative_x) .* (target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)])) ) ;x_scaling(2) =(target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)])) ) ;

y_scaling = (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)]))
y_rotated_scaled = y_rotated .* y_scaling;

%% 2: thes eare the dots from the GUI just for orientation
maze.NE_rot_scaled(1) = maze.NE_rot(1) .* target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)]));maze.SE_rot_scaled(1) = maze.SE_rot(1) .* target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)])); 
maze.NW_rot_scaled(1) = -maze.NW_rot(1) .* target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)]));maze.SW_rot_scaled(1) = -maze.SW_rot(1) .* target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)]));
maze.NE_rot_scaled(2) = maze.NE_rot(2)* (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)])); maze.NW_rot_scaled(2) = maze.NW_rot(2) * (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)]));
maze.SE_rot_scaled(2) = maze.SE_rot(2)* (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)])); maze.SW_rot_scaled(2) = maze.SW_rot(2) * (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)]));


figure
scatter([maze.SW_rot_scaled(1), maze.SE_rot_scaled(1), maze.NW_rot_scaled(1), maze.NE_rot_scaled(1)], [maze.SW_rot_scaled(2), maze.SE_rot_scaled(2), maze.NW_rot_scaled(2), maze.NE_rot_scaled(2)], 'g')
hold on
plot(x_rotated_scaled, y_rotated_scaled, 'r-', 0, 0, 'bo');
axis equal
hold off
title('Rotated and scaled');
saveas(gcf,fullfile(save_dir,'Rotated_AndScaled.png')); 
close

%-----------------------------------------------------
%% Scaling for Sleep
%-----------------------------------------------------

aver_x_scaling = mean(abs(x_scaling));
scaling_for_sleep = mean(abs([aver_x_scaling, y_scaling]));

%% Rotate 180 degrees for linearization of the path
%{
theta = deg2rad(180);
v = [x_rotated_scaled; y_rotated_scaled];

% choose a point which will be the center of rotation

x_center = 0;

y_center = 0;

% create a matrix which will be used later in calculations

center = repmat([x_center; y_center], 1, length(x));
%angle = 
%rotation_matrix = repmat([

% define a 60 degree counter-clockwise rotation matrix

%theta = pi/3;       % pi/3 radians = 60 degrees

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% do the rotation...

s = v - center;     % shift points in the plane so that the center of rotation is at the origin

so = R*s;  

x_final = so(1,:);

y_final = so(2,:);
figure
plot(x_final,y_final);
title('180 degrees rotated');
%}

end

%% 3: Scaling arms increasingly, not a good idea
%{
x_rotated_test2(positive_x) = (x_rotated(positive_x)).^2 .* (26 / (mean([maze.NE_rot(1) maze.SE_rot(1)])).^2 ) ; 
x_rotated_test2(negative_x) = -(x_rotated(negative_x)).^2 .* (26 / (mean([maze.NW_rot(1) maze.SW_rot(1)])).^2 ) ;
figure

plot(x_rotated_test2, y_rotated, 'r-', 0, 0, 'bo');
hold off
%}

