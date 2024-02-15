function [locseq, succArms] = fixsuccarmseq(locseq)
% put a node in between two successive arms

locseqin = locseq;

succArms = cellfun(@(lag,lead) ~isempty(regexp(lag, 'A', 'once')) && strcmp(regexp(lag, 'A', 'match', 'once'), regexp(lead, 'A', 'match', 'once')), locseq(1:end-1), locseq(2:end));
succArms = [false; succArms];

fixedseq = repmat({''}, length(locseq)+sum(succArms), 1);
holes = find(succArms);
succArms = num2logical(holes + find(holes) - 1, length(fixedseq));
fixedseq(~succArms) = locseq;
fixedseq = fillSeqHole(fixedseq);

while length(fixedseq) ~= length(locseq)
	locseq = fixedseq;
	fixedseq = fixsuccarmseq(locseq);
end
locseq = fixedseq;

succArms = ones(length(locseq), 1);
j = 1;
for i = 1:length(locseq)
	if strcmp(locseqin{j}, locseq{i})
		succArms(i) = 0;
		j = j + 1;
	end
end