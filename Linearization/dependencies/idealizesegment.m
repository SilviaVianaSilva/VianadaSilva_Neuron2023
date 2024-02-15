function [x y] = idealizesegment(x0,y0,xNode,yNode,entryExit)

% comments added by @author: Sia @date 9/15/15 5:27 PM

if all(~ismember({'l','r'},entryExit)) % a side arm @author: Sia @date 9/15/15
    x = x0;
    y = y0;
    x(:) = mean(minmax(xNode)); % gets rid of horizontal displacement
elseif all(~ismember({'t','b'},entryExit)) % a base arm @author: Sia @date 9/15/15
    x = x0;
    y = y0;
    y(:) = mean(minmax(yNode)); % gets rid of vertical displacement
else
    if ismember('t',entryExit)
        dy = max(yNode);
    else
        dy = min(yNode);
    end
    if ismember('l',entryExit)
        dx = min(xNode);
    else
        dx = max(xNode);
    end
    a = cart2pol(x0-dx,y0-dy);
    r = diff(minmax(xNode))/2; % half length middle arm
    [x y] = pol2cart(a,r); % Is he scaling?
    x = x+dx;
    y = y+dy;
end
        