function [rect_x, rect_y] = rectanch(w,h,anchor,scale)

if ~exist('scale', 'var') || isempty(scale)
	scale = 1;
end

if ~exist('anchor', 'var') || isempty(anchor)
	anchor = 1;
end



rect_x = [0 0 w w 0] * scale;
rect_y = [0 h h 0 0] * scale;

switch anchor
	case 2
		rect_y = rect_y - h;
	case 3
		rect_x = rect_x - w;
		rect_y = rect_y - h;
	case 4
		rect_x = rect_x - w;
	case 0
		rect_x = rect_x - w/2;
		rect_y = rect_y - h/2;
	otherwise
		% pass
end

rect_x = rect_x(:);
rect_y = rect_y(:);