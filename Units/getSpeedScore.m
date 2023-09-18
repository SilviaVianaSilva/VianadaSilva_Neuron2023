function [speedScore, pVal] = getSpeedScore(time,velocity,tSp,minSpeed)
if ~exist('minSpeed','var') || isempty(minSpeed)
    minSpeed = 2;
end
[~, IFR] = getInstFR(time,tSp);
velocity(velocity<minSpeed) = nan;
if length(tSp)<50
    speedScore = nan;
    pVal = nan;
else
    [speedScore, pVal] = nancorr_speedmod(IFR,velocity');
end






