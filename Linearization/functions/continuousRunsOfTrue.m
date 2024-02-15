function beginAndEnd = continuousRunsOfTrue(vectorOfLogicals)
% SYNTAX: beginAndEnd = continuousRunsOfTrue(vectorOfLogicals)
% returns an n-by-2 array of indices, such that for each row, the first
% element is the start of an uninterrupted run of true and the second 
% element is the end.
% Note that a single true value at position j will return a row [j j].
% If there are no true values in vectorOfLogicals, this returns [0 0]
%
% Emily Mankin

% convert to ones and zeros, just in case...
vectorOfLogicals = double(logical(vectorOfLogicals));

% special cases:
if sum(vectorOfLogicals)==0
    beginAndEnd = zeros(0,2);
elseif sum(vectorOfLogicals)==length(vectorOfLogicals)
    beginAndEnd = [1,length(vectorOfLogicals)];
else
    % normal cases:
    v = vectorOfLogicals;
    dv = diff(v);
    endPoints = find(dv==-1);
    if v(end)==1&&(isempty(endPoints)||endPoints(end)~=length(v));
        endPoints = [endPoints,length(v)];
    end
    startPoints = find(dv==1)+1;
    if isempty(startPoints) || endPoints(1)<startPoints(1)
        startPoints = [1,startPoints];
    end
    if isempty(endPoints)||startPoints(end)>endPoints(end)
        endPoints = [endPoints,length(v)];
    end
    beginAndEnd = [startPoints',endPoints'];
end