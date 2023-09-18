function [x,y,t,angle,n_samples, n_trackingPoints,n_NaNPoints] = getPreprocessedPositionData(videoFile,doublediode)
if doublediode == 1
    [t, x, y, angle, n_samples] = writeVelocity.ProcessVideoData_double(videoFile,videoFile(1:end-7));%
    n_trackingPoints = n_samples;
    n_NaNPoints = length(find(isnan(x)));
else
    [t, x, y, angle, n_samples] = writeVelocity.ProcessVideoData(videoFile,doublediode);
    n_trackingPoints = n_samples;
    n_NaNPoints = n_samples - length(x);
end
index = find(~isnan(x));
t = t(index);
x = x(index);
y = y(index);
angle = angle(index);

t = t/1000000;

