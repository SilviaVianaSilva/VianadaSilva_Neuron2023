function d = getPathDistance(x,y,type,x0,y0)
if nargin<3 || isempty(type)
    type = 'trajectory';
end

switch type
    case 'trajectory'
        dx = diff(x);
        dy = diff(y);
        dd = zeros(length(x),1);
        for k = 1:length(dx)
            dd(k+1) = norm([dx(k) dy(k)]);
        end
        d = cumsum(dd);
    case 'absolute'
        if ~exist('x0','var')
            x0 = x(1);
        end
        if ~exist('y0','var')
            y0 = y(1);
        end            
        dx = x-x0;
        dy = y-y0;
        d = zeros(length(x),1);
        for k = 1:length(dx)
            d(k) = norm([dx(k) dy(k)]);
        end
end