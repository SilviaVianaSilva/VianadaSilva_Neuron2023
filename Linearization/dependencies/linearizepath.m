function pathDataLin = linearizepath(pathDataIdeal,trialInfo)

nPaths = length(pathDataIdeal);
pathDataLin = pathDataIdeal;
for path = 1:nPaths
    mazeType = trialInfo(path,1).mazeType;
    switch mazeType
        case {'fig8', 'fig8rat', 'fig8mouse'}
        template = fig8trialtemplate2(mazeType);
    end
    nTrials = size(trialInfo(path).inds,1);
    x = pathDataIdeal(path).x;
    y = pathDataIdeal(path).y;
    xLin = nan(size(x));
    yLin = nan(size(y));
    for tr = 1:nTrials
       dir = trialInfo(path).direction{tr}; 
       inds = trialInfo(path).inds{tr}(1):trialInfo(path).inds{tr}(end);
       xLin(inds) = linterp2(template.(dir).x,template.(dir).y,template.(dir).d,x(inds),y(inds));
       yLin(inds) = 0;
    end
    pathDataLin(path,1).x = xLin;
    pathDataLin(path,1).y = yLin;
end