function [locseq, succNodes] = fixsuccnodeseq(locseq)
% put the arm in between two successive arms

locseqin = locseq;

succNodes = cellfun(@(lag,lead) ~isempty(regexp(lag, 'N', 'once')) && strcmp(regexp(lag, 'N', 'match', 'once'), regexp(lead, 'N', 'match', 'once')), locseq(1:end-1), locseq(2:end));
succNodes = [false; succNodes];

fixedseq = repmat({''}, length(locseq)+sum(succNodes), 1);
holes = find(succNodes);
succNodes = num2logical(holes + find(holes) - 1, length(fixedseq));
fixedseq(~succNodes) = locseq;
fixedseq = fillArmHole(fixedseq);

while length(fixedseq) ~= length(locseq)
	locseq = fixedseq;
	fixedseq = fixsuccarmseq(locseq);
end
locseq = fixedseq;

succNodes = ones(length(locseq), 1);
j = 1;
for i = 1:length(locseq)
	if strcmp(locseqin{j}, locseq{i})
		succNodes(i) = 0;
		j = j + 1;
	end
end