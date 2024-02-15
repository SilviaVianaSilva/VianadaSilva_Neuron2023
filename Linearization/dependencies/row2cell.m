function C = row2cell(S)

if ~ismatrix(S)
	error('S must be a 2-d array')
end

C = mat2cell(S, ones(size(S, 1), 1), size(S, 2));