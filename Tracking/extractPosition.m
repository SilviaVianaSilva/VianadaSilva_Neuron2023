% Extract the x and y coordinates from the target data, while applying
% constrainst on the target. diodeColor defines the color to be extraced.
function [x,y] = extractPosition(targets,diodeColor)
if ischar(diodeColor)
    diodeColor = {diodeColor};
end
nColors = length(diodeColor);
% Number of samples in the position file
nSamples = size(targets,1);
x = zeros(nSamples,nColors);
y = zeros(nSamples,nColors);
for c = 1:nColors
    currentTargets = targets;
    % Take out targets for colours not in use
    switch diodeColor{c}
        case {'l','lum','luminance'}
            % Keep only luminance
            for ii = 1:nSamples
                index = ~currentTargets(ii,:,3);
                currentTargets(ii,index,:) = 0;
            end
        case {'r','red'}
            % Keep only red
            for ii = 1:nSamples
                index = ~currentTargets(ii,:,4) & ~currentTargets(ii,:,7);
                currentTargets(ii,index,:) = 0;
            end
        case {'g','green'}
            % Keep only green
            for ii = 1:nSamples
                index = ~currentTargets(ii,:,5) & ~currentTargets(ii,:,8);
                currentTargets(ii,index,:) = 0;
            end
        case {'b','blue'}
            % Keep only blue
            for ii = 1:nSamples
                index = ~currentTargets(ii,:,6) & ~currentTargets(ii,:,9);
                currentTargets(ii,index,:) = 0;
            end
    end
    
    % Calculate the position for the first sample(s)
    sx = [];
    sy = [];
    currentSamp = 1;
    while currentSamp <= nSamples
        % Index to the legal targets for this sample
        tarInd = find(currentTargets(currentSamp,:,1)>0);
        % Count number of targets for the current sample
        n = length(tarInd);
        
        if n > 0
            % Set the start position to the mean of the targets
            sx = mean(currentTargets(currentSamp,tarInd,1));
            sy = mean(currentTargets(currentSamp,tarInd,2));
            % First sample position is found and the loop is terminated
            break
        end
        % Set index for the next sample candidate
        currentSamp = currentSamp + 1;
    end
    
    % No targets for this color
    if isempty(sx)
        return
    end
    
    x(currentSamp,c) = sx;
    y(currentSamp,c) = sy;
    
    % Go through the rest of the targets and extract the positions
    for ii = currentSamp+1:nSamples
        % Index to the legal targets for this sample
        tarInd = find(currentTargets(ii,:,1)>0);
        % Find number of currentTargets for this sample (max 50)
        n = length(tarInd);
        if n > 0
            % Sample position is set to the mean of the legal targets
            x(ii,c) = mean(targets(ii,tarInd,1));
            y(ii,c) = mean(targets(ii,tarInd,2));
        end
    end
end
