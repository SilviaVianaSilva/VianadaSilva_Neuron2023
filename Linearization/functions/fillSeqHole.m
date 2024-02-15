function lblSeq = fillSeqHole(lblSeq)

% 03/02/2016 8:59 PM
% Siavash Ahmadi

correctSeq = {'A34', 'A23', 'A25', 'A56', 'A16', 'A12', 'A25', 'A45'};


idx = strcmpi(lblSeq, '');
idx(1) = false; idx(end) = false; % for the first and last indices we can never be sure what the previous location of the animal has been (unless we look at the videoData)

for i = find(idx(:)')
	arm1 = lblSeq{i-1};
	arm2 = lblSeq{i+1};
	buffer = setdiff(intersect(arm1, arm2), 'A'); % the common number in the two successive elements
	if isempty(buffer) % non-adjacent arms
		possibleArms = armsbetween(arm1, arm2);
		if isempty(possibleArms) % 2-adjacent arms
			% TODO
			buffer = '';
		else % 1-adjacent arms
			[~, ~, ic] = intersect(possibleArms, correctSeq);
			% the following line chooses the arm that will ensure a wrong sequence
			% (this is so that the trial will be marked degenerate because it is
			% ambiguous which arm is the true arm)
			buffer = possibleArms(strcmp(correctSeq(mod(ic, length(correctSeq))+1), arm1));
		end
		lblSeq(i) = buffer;
	else
		lblSeq{i} = ['N', buffer];
	end
end

function possibleArms = armsbetween(arm1, arm2)
allArms = {'A16';
	'A12';
	'A25';
	'A23';
	'A45';
	'A34';
    'A56'};
arm1nodes = finite(str2double(num2cell(arm1)));
arm2nodes = finite(str2double(num2cell(arm2)));
possibleArms = [arm1nodes(1), arm2nodes(1);
	arm1nodes(1), arm2nodes(2);
	arm1nodes(2), arm2nodes(1);
	arm1nodes(2), arm2nodes(2);];
possibleArms = [min(possibleArms, [], 2), max(possibleArms, [], 2)];
possibleArms = strcat('A', num2str(possibleArms, '%d%d'));
possibleArms = row2cell(possibleArms);
possibleArms = intersect(possibleArms, allArms);