function [IFR, smoothedIFR] = getInstFR(time,tSp,smoothingWindow)

if ~exist('smoothingWindow','var') || isempty(smoothingWindow)
    smoothingWindow = .25;
end

points = length(time);
IFR = nan(points-1,1);

for i = 1:points-1
    start = time(i);
    finish = time(i+1);
    nSpikes = length(find(iswithin(tSp,start,finish)));
    duration = finish-start;
    IFR(i) = nSpikes/duration;
end

meanFrameSize = median(diff(time));
smoothBinWindow = round(smoothingWindow/meanFrameSize);

gaussian = gausswin(smoothBinWindow);

smoothedIFR = conv(IFR,gaussian);

IFR = IFR((smoothBinWindow+1)/2:end-(smoothBinWindow-1)/2);
smoothedIFR = smoothedIFR((smoothBinWindow+1)/2:end-(smoothBinWindow-1)/2);

IFR = [nan IFR'];
smoothedIFR = [nan smoothedIFR'];