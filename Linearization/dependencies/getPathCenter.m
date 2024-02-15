function [xC yC] = getPathCenter(x,y)
maxX = max(x);
minX = min(x);
maxY = max(y);
minY = min(y);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
center = findCentre(NE,NW,SW,SE);
xC = center(1);
yC = center(2);