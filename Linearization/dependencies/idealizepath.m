function [pathDataIdeal, locInfo, maze] = idealizepath(pathData,mazeType)
if nargin<2
    mazeType = 'fig8';
end
[locInfo, maze] = mazeLocInfo(pathData,mazeType);
nPaths = length(pathData);
pathDataIdeal = pathData;
for path = 1:nPaths
    labelSeq = locInfo(path).labelSeq;
    locInds = locInfo(path).inds;
	locInfo(path).mazeType = mazeType;
    nSeq = length(labelSeq);
    x = pathData(path).x;
    y = pathData(path).y;
    for i = 1:nSeq
        loc = labelSeq{i};
        inds = locInds(i,1):locInds(i,2);
        if isempty(loc)
            x(inds) = nan;
            y(inds) = nan;
            continue;
        end
        if i == 1
            entryExitLoc = {labelSeq{i+1},labelSeq{i+1}};
        elseif i == nSeq
            entryExitLoc = {labelSeq{i-1},labelSeq{i-1}};
        else
            entryExitLoc = {labelSeq{i-1},labelSeq{i+1}};
        end
        xBnd = maze.locs.(loc)(:,1);
        yBnd = maze.locs.(loc)(:,2);
        switch loc
            case {'A16','A25','A34'}
                entryExit = {'t','b'};
            case {'A12','A23','A45','A56'};
                entryExit = {'l','r'};
            case 'N1'
                entryExit = {'b','r'};
            case 'N2'
                code = {'l','b','r'};
                adj = {'A12','A25','A23'};
                entryExit = code(ismember(adj,entryExitLoc));
            case 'N3'
                entryExit = {'l','b'};
            case 'N4'
                entryExit = {'t','l'};

            case 'N5'
                code = {'r','t','l'};
                adj = {'A45','A25','A56'};
                entryExit = code(ismember(adj,entryExitLoc));
            case 'N6'
                entryExit = {'r','t'};
        end
        [x(inds), y(inds)] = idealizesegment(x(inds),y(inds),xBnd,yBnd,entryExit);
    end
    pathDataIdeal(path,1).x = x;
    pathDataIdeal(path,1).y = y;
end