function l = num2logical(numbers, mxlen)

if nargin < 2
	mxlen = max(numbers(:));
end

if isvector(numbers)
	if isrow(numbers)
		l = false(1, mxlen);
	else
		l = false(mxlen, 1);
	end
else
	l = false(mxlen, 1);
end

l(numbers(:)) = true;