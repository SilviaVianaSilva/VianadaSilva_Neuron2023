function out = output(N, f, varargin)
%OUTPUT Capture N-th output of @f with inputs X1,...,Xn where 1 <= N <= n
%
% out = OUTPUT(N, f, varargin)

% Siavash Ahmadi
% 3/4/2016 11:15 AM

[x{1:nargout(f)}] = f(varargin{:});
out = x{N};