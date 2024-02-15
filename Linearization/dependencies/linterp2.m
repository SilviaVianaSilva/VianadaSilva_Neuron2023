function zi = linterp2(x,y,z,xi,yi)
warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi);
warning on MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
% [z iUnq] = unique(z);
% x = x(iUnq);
% y = y(iUnq);
% n = length(xi);
% zi = nan(n,1);
% dzMax = max(diff(z));
% warning off;
% for i = 1:n
%     R =  arrayfun(@(j)norm([xi(i)-x(j) yi(i)-y(j)]),1:length(x));
%     [~, iSort] = sort(R);
%     c = sort(iSort(1:2));
%     if all(c==c(1)) || diff(c)~=1 % ||max(R(c))>dzMax
%         zi(i) = nan;
%     else
%         zi(i) = mean([z(c(1))+R(c(1)) z(c(2))-R(c(2))]);
% %         m = [x(c) y(c)]\z(c);
% %         zi(i) = xi(i)*m(1) + yi(i)*m(2);
%         if ~in(zi(i),z(c))
%             zi(i) = nan;
%         end
%     end
% end
% warning on;