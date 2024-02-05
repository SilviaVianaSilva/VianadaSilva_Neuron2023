%Maze_Preview_MH
target_maze = [52 72]; % Insert here the total length of the maze

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


fh2 = figure()
scatter([maze.SW_rot_scaled(1), maze.SE_rot_scaled(1), maze.NW_rot_scaled(1), maze.NE_rot_scaled(1)], [maze.SW_rot_scaled(2), maze.SE_rot_scaled(2), maze.NW_rot_scaled(2), maze.NE_rot_scaled(2)], 'g')
hold on
plot(x_rotated_scaled, y_rotated_scaled, 'r-', 0, 0, 'bo');
axis equal
hold off
title('Rotated and scaled');

