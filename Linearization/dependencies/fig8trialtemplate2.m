function [template, locBnds] = fig8trialtemplate2(mazeType)
% Modified from Chris's fig8trialtemplate to produce a different trial
% progression (starting at post-reward and ending at reward)

% 3/3/2016
% Siavash Ahmadi

if nargin == 0
	mazeType = 'fig8rat';
end

maze = fig8maze(mazeType);
locs = fields(maze.locs);
[xC, yC] = getPathCenter(maze.whole.basic(:,1),maze.whole.basic(:,2));
for l = 1:length(locs)
   maze.locs.(locs{l})(:,1) = maze.locs.(locs{l})(:,1)-xC;
   maze.locs.(locs{l})(:,2) = maze.locs.(locs{l})(:,2)-yC;
end
r = diff(minmax(maze.locs.N1(:,1)))/2;
dd = 0.2;
da = dd/r;


xTmp = [];
yTmp = [];
dTmp = 0;
lastMax = 0;


yA16 = [min(maze.locs.A16(:,2)):dd:max(maze.locs.A16(:,2))]';
xA16 = sum(minmax(maze.locs.A16(:,1)))/2*ones(length(yA16),1);

[xN1, yN1] = pol2cart(-[pi:da:3*pi/2]',r);
xN1 = xN1 + max(maze.locs.N1(:,1));
yN1 = yN1 + min(maze.locs.N1(:,2));

xA12 = [min(maze.locs.A12(:,1)):dd:max(maze.locs.A12(:,1))]';
yA12 = sum(minmax(maze.locs.A12(:,2)))/2*ones(length(xA12),1);

xTmp = cat(1,xTmp,[xA16;xN1;xA12]);
yTmp = cat(1,yTmp,[yA16;yN1;yA12]);
dTmp = cat(1, [], max(dTmp) + getPathDistance([xA16;xN1;xA12],[yA16;yN1;yA12]), max(dTmp) + getPathDistance([xA16;xN1;xA12],[yA16;yN1;yA12]));
xTmp = cat(1,xTmp,-xTmp);
yTmp = cat(1,yTmp,yTmp);

locBnds.return = minmax(dTmp(dTmp>=lastMax));
lastMax = max(dTmp);




[xN2L, yN2L] = pol2cart(-[3*pi/2:da:2*pi]',r);
xN2L = xN2L + min(maze.locs.N2(:,1));
yN2L = yN2L + min(maze.locs.N2(:,2));

xTmp = cat(1,xTmp,xN2L);
yTmp = cat(1,yTmp,yN2L);
xTmp = cat(1,xTmp,-xN2L);
yTmp = cat(1,yTmp,yN2L);
dTmp = cat(1, dTmp, max(dTmp) + getPathDistance(xN2L,yN2L), max(dTmp) + getPathDistance(xN2L,yN2L));

locBnds.delay = minmax(dTmp(dTmp>=lastMax));
lastMax = max(dTmp);



yA25 = sort([min(maze.locs.A25(:,2)):dd:max(maze.locs.A25(:,2))]','descend');
xA25 = sum(minmax(maze.locs.A25(:,1)))/2*ones(length(yA25),1);

xTmp = cat(1,xTmp,xA25);
yTmp = cat(1,yTmp,yA25);
dTmp = cat(1, dTmp, max(dTmp) + getPathDistance(xA25, yA25));
locBnds.center = minmax(dTmp(dTmp>=lastMax));
lastMax = max(dTmp);

[xN5L, yN5L] = pol2cart(-[0:da:pi/2]',r);
xN5L = xN5L + min(maze.locs.N5(:,1));
yN5L = yN5L + max(maze.locs.N5(:,2));
xTmp = cat(1,xTmp,xN5L);
yTmp = cat(1,yTmp,yN5L);
dTmp = cat(1, dTmp, max(dTmp) + getPathDistance(xN5L, yN5L));
locBnds.bifur = minmax(dTmp(dTmp>=lastMax));
lastMax = max(dTmp);

[xN5R, yN5R] = pol2cart([pi:da:5*pi/2]',r);
xN5R = xN5R + max(maze.locs.N5(:,1));
yN5R = yN5R + max(maze.locs.N5(:,2));

xA56 = sort([min(maze.locs.A56(:,1)):dd:max(maze.locs.A56(:,1))]','descend');
yA56 = sum(minmax(maze.locs.A56(:,2)))/2*ones(length(xA56),1);
xTmp = cat(1,xTmp,xA56);
yTmp = cat(1,yTmp,yA56);
dTmp = cat(1, dTmp, max(dTmp) + getPathDistance(xA56, yA56));
locBnds.rewardarm = minmax(dTmp(dTmp>=lastMax));
lastMax = max(dTmp);

[xN6, yN6] = pol2cart(-[pi/2:da:pi]',r);
xN6 = xN6 + max(maze.locs.N6(:,1));
yN6 = yN6 + max(maze.locs.N6(:,2));
xTmp = cat(1,xTmp,xN6);
yTmp = cat(1,yTmp,yN6);
dTmp = cat(1, dTmp, max(dTmp) + getPathDistance(xN6, yN6));
locBnds.reward = minmax(dTmp(dTmp>=lastMax));
lastMax = max(dTmp);







d = dTmp;

template.L.x = xTmp;
template.L.y = yTmp;
template.L.d = d;
template.R.x = -xTmp;
template.R.y = yTmp;
template.R.d = d;
    
