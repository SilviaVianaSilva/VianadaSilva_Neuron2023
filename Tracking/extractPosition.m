% Extract the x and y coordinates from the target data, while applying
% constrainst on the target. tarCol defines the colour to be extraced.
function [x,y,targets,exFlag] = extractPosition(targets,borderX,borderY,tarCol)

% Convert border values to int16 (same as targets)
borderX = int16(borderX);
borderY = int16(borderY);

% Number of samples in the position file
numSamp = size(targets,1);
x = zeros(numSamp,1);
y = zeros(numSamp,1);

% Take out targets for colours not in use
switch tarCol
    case 1
        % Keep only luminance
        for ii = 1:numSamp
            index = find(~targets(ii,:,3));
            targets(ii,index,:) = 0;
        end
    case 2
        % Keep only red
        for ii = 1:numSamp
            index = find(~targets(ii,:,4) & ~targets(ii,:,7));
            targets(ii,index,:) = 0;
        end
    case 3
        % Keep only green
        for ii = 1:numSamp
            index = find(~targets(ii,:,5) & ~targets(ii,:,8));
            targets(ii,index,:) = 0;
        end
    case 4
        % Keep only blue
        for ii = 1:numSamp
            index = find(~targets(ii,:,6) & ~targets(ii,:,9));
            targets(ii,index,:) = 0;
        end
end

% % Remove every target that are outside the border set by the user
% for ii = 1:numSamp
%     inside = inpolygon(targets(ii,:,1),targets(ii,:,2),borderX,borderY);
%     targets(ii,inside==0,1:2) = 0;
% end

% Calculate the position for the first sample(s)
sx = [];
sy = [];
currentSamp = 1;
while currentSamp <= numSamp
    % Index to the legal targets for this sample
    tarInd = find(targets(currentSamp,:,1)>0);
    % Count number of targets for the current sample
    n = length(tarInd);

    if n > 0
        % Set the start position to the mean of the targets
        sx = mean(targets(currentSamp,tarInd,1));
        sy = mean(targets(currentSamp,tarInd,2));
        % First sample position is found and the loop is terminated
        break
    end
    % Set index for the next sample candidate
    currentSamp = currentSamp + 1;
end

% No targets for this colour
if isempty(sx)
    exFlag = 0;
    return
end

x(currentSamp) = sx;
y(currentSamp) = sy;

% Go through the rest of the targets and extract the positions
for ii = currentSamp+1:numSamp
    % Index to the legal targets for this sample
    tarInd = find(targets(ii,:,1)>0);
    % Find number of targets for this sample (max 50)
    n = length(tarInd);

    if n > 0
        % Sample position is set to the mean of the legal targets
        x(ii) = mean(targets(ii,tarInd,1));
        y(ii) = mean(targets(ii,tarInd,2));
    end
end

exFlag = 1;
