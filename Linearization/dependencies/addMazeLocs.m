function [pathData, maze] = addMazeLocs(pathData,mazeType)
if nargin<2
    mazeType = 'fig8mouse';
end
switch mazeType
    case {'fig8', 'fig8rat', 'fig8mouse'}
        maze = fig8maze(mazeType);
        
end
[xC, yC] = getPathCenter(maze.boundary.outer(:,1),maze.boundary.outer(:,2));
mazeFields = fields(maze);
for i = 1:length(mazeFields)
    subFields = fields(maze.(mazeFields{i}));
    for j = 1:length(subFields)
        maze.(mazeFields{i}).(subFields{j})(:,1) = maze.(mazeFields{i}).(subFields{j})(:,1)-xC;
        maze.(mazeFields{i}).(subFields{j})(:,2) = maze.(mazeFields{i}).(subFields{j})(:,2)-yC;
    end
end
nPaths = length(pathData);
for path = 1:nPaths
    nPts = length(pathData(path).x);
    locID = zeros(nPts,1);
    locs = fields(maze.locs);
    nLocs = length(locs);
    
    
  
    for i = 1:nLocs
        locInds = inpolygon(pathData(path).x,pathData(path).y,maze.regx.(locs{i})(:,1), maze.regx.(locs{i})(:,2));
        locID(locInds) = i;
    end
	
     
	%%% Sia's fixing algorithm
	locID(1) = locID(find(locID, 1)); % let the first element of locID be the first non-zero element
	for i = 2:length(locID)-1
		if locID(i) == 0
			if locID(i-1) == locID(i+1)
				locID(i) = locID(i-1);
			elseif locID(i+1) > 6
				locID(i) = str2double(intersect(intersect(locs{locID(i-1)}, locs{locID(i+1)}), num2str(1:6)));
			else
				locID(i) = locID(i-1);
			end
		elseif locID(i) > 6 % an arm
			if locID(i+1) > 6 && locID(i+1) ~= locID(i)
				locID(i) = str2double(intersect(intersect(locs{locID(i)}, locs{locID(i+1)}), num2str(1:6)));
			end
		else % a node
			if locID(i+1) <= 6 && locID(i+1) >= 1 && locID(i) ~= locID(i+1) % node-to-node transition
				m = min(locID(i), locID(i+1));
				M = max(locID(i), locID(i+1));
				locID(i) = find(strcmp(locs, ['A', num2str(m), num2str(M)]));
			end
		end
	end
	locID(isnan(locID)) = 0;
	%%%
	
    noLoc = locID==0;
    noLocBnds = continuousRunsOfTrue(noLoc');
    nBnds = size(noLocBnds,1);
    for b = 1:nBnds
        prevLoc = locID(max(1,noLocBnds(b,1)-1));
		if noLocBnds(b,2) < length(noLocBnds(b,2))
			nextLoc = locID(noLocBnds(b,2)+1);
		else
			nextLoc = prevLoc;
		end
        if prevLoc == nextLoc
            locID(noLocBnds(b,1):noLocBnds(b,2)) = prevLoc;
        end
    end
    pathData(path).locID = locID;
    pathData(path).locLabel = cell(nPts,1);
    hasLoc = locID>0;
    pathData(path).locLabel(hasLoc) = locs(locID(hasLoc));
    pathData(path).locLabel(~hasLoc) = {''};
end