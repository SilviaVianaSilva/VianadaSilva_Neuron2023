function lblSeq = fillArmHole(lblSeq)

% 10/21/2016 3:57 PM
% Siavash Ahmadi

correctSeqArm = {'A34', 'A23', 'A25', 'A56', 'A16', 'A12', 'A25', 'A45'};
correctSeqNodes = {'N3', 'N2', 'N5', 'N6', 'N1', 'N2', 'N5', 'N4'};

idx = strcmpi(lblSeq, '');
idx(1) = false; idx(end) = false; % for the first and last indices we can never be sure what the previous location of the animal has been (unless we look at the videoData)

for i = find(idx(:)')
	node1 = lblSeq{i-1};
	node2 = lblSeq{i+1};
	buffer = sort([setdiff(node1, 'N'), setdiff(node2, 'N')]); % the common number in the two successive elements
	if length(buffer) ~= 2 % non-adjacent arms
		possibleArms = armsbetween(node1, node2);
		if isempty(possibleArms) % 2-adjacent arms
			% TODO
			buffer = '';
		else % 1-adjacent arms
			[~, ~, ic] = intersect(possibleArms, correctSeq);
			% the following line chooses the arm that will ensure a wrong sequence
			% (this is so that the trial will be marked degenerate because it is
			% ambiguous which arm is the true arm)
			buffer = possibleArms(strcmp(correctSeq(mod(ic, length(correctSeq))+1), node1));
		end
		lblSeq(i) = buffer;
	else
		lblSeq{i} = ['A', buffer];
	end
end