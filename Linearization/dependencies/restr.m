function [r, I, iNaN] = restr(S, low, hi, sepByNaN)
%[r, I, iNaN] = restr(S, low, hi, sepByNaN)

% Siavash Ahmadi
% 03/22/2015 Version 1.0
% 12/13/2015 Version 2.0 (added functionality: handles discontguous
% intervals, uses NaNs to demarcate intervals, preserves input-output shape)

if ~exist('sepByNaN', 'var') || isempty(sepByNaN)
	sepByNaN = false;
end

iNaN = [];
withinEachRange = arrayfun(@(l,h) l <= S & S <= h, low, hi, 'un', 0);
I = false(size(S));
for i = 1:length(withinEachRange)
	I = I | withinEachRange{i};
end

if sepByNaN
	if sum(cellfun(@any, withinEachRange)) == 1 % don't put a NaN at the end of a contiguous array
		r = S(withinEachRange{cellfun(@any, withinEachRange)});
		iNaN = false(size(r));
		return;
	end
	r_cell = cellfun(@(idx) accFunc_ExtractAndAddNaN(S,idx), withinEachRange, 'un', 0);
	r = cat(1, r_cell{:});
	iNaN = false(size(r));
	iNaN(cumsum(cellfun(@length, r_cell))) = true;
	% TODO: make sure shape of |r| is consistent with that of S
else
	r = S(I);
end

function S = accFunc_ExtractAndAddNaN(S,idx)
S = S(idx);
S(end+1) = NaN;
S = S(:);