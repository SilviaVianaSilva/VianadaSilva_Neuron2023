function [x y] = rectcoords(w,h,x0,y0)
if nargin<4
    y0 = 0;
    x0 = 0;
end
if nargin<2
    h = w;
end

x = w*[0 1 1 0 0]' + x0;
y = h*[1 1 0 0 1]' + y0;