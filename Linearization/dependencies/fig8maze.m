function maze = fig8maze(w,h,wNode,hNode)
CM = 2.54;

if nargin == 0
    mazeType = 'fig8rat';
end
if nargin == 1
	mazeType = w;
end

if strcmpi(mazeType, 'fig8rat') || strcmpi(mazeType, 'fig8')
	w = (40+5/8)*CM;
	h = (5*12)*CM;
elseif strcmpi(mazeType, 'fig8mouse')
	w = 20*CM;
	h = 30*CM;
elseif ischar(mazeType)
    error('Bad mazeType');
end

if nargin<4
    wNode = w*5/(40+5/8);
    hNode = h/12;
end
lenHArm = (w-3*wNode)/2;
lenVArm = h-2*hNode;

[xOuter, yOuter] = rectcoords(w,h);
[xInner, yInner] = rectcoords(w-2*wNode,h-2*wNode);
xInner = xInner + wNode;
yInner = yInner + hNode;

maze.boundary.outer = [xOuter yOuter];
maze.boundary.inner = [xInner yInner];

maze.whole.basic = [xOuter yOuter; nan nan; xInner yInner];

[xNode, yNode] = rectcoords(wNode,hNode);
[xVArm, yVArm] = rectcoords(wNode,lenVArm);
[xHArm, yHArm] = rectcoords(lenHArm,hNode);

maze.locs.N1 = [xNode yNode+h-hNode];
maze.locs.N2 = [xNode+wNode+lenHArm yNode+h-hNode];
maze.locs.N3 = [xNode+w-wNode yNode+h-hNode];
maze.locs.N4 = [xNode+w-wNode yNode];
maze.locs.N5 = [xNode+wNode+lenHArm yNode];
maze.locs.N6 = [xNode yNode];


stickoutMultiplier = 3;

[xMain, yMain] = rectanch(wNode*stickoutMultiplier,hNode*stickoutMultiplier,4);
[xAux, yAux] = rectanch(lenHArm/2,hNode,2);
maze.regx.N1 = [xMain+wNode, yMain+h-hNode];
maze.regx.N1 = [maze.regx.N1; [NaN NaN]; xAux+wNode, yAux+h-hNode];

[xMain, yMain] = rectanch(wNode,hNode*stickoutMultiplier,1);
[xAux, yAux] = rectanch(lenHArm/2,hNode,3);
maze.regx.N2 = [xMain+wNode+lenHArm, yMain+h-hNode];
maze.regx.N2 = [maze.regx.N2; [NaN NaN]; xAux+wNode+lenHArm, yAux+h-hNode];
[xAux, yAux] = rectanch(lenHArm/2,hNode,2);
maze.regx.N2 = [maze.regx.N2; [NaN NaN]; xAux+wNode*2+lenHArm, yAux+h-hNode];

[xMain, yMain] = rectanch(wNode*stickoutMultiplier,hNode*stickoutMultiplier,1);
[xAux, yAux] = rectanch(lenHArm/2,hNode,3);
maze.regx.N3 = [xMain+wNode*2+lenHArm*2, yMain+h-hNode];
maze.regx.N3 = [maze.regx.N3; [NaN NaN]; xAux+wNode*2+lenHArm*2, yAux+h-hNode];

[xMain, yMain] = rectanch(wNode*stickoutMultiplier,hNode*stickoutMultiplier,2);
[xAux, yAux] = rectanch(lenHArm/2,hNode,4);
maze.regx.N4 = [xMain+wNode*2+lenHArm*2, yMain+hNode];
maze.regx.N4 = [maze.regx.N4; [NaN NaN]; xAux+wNode*2+lenHArm*2, yAux+hNode];

[xMain, yMain] = rectanch(wNode,hNode*stickoutMultiplier,2);
[xAux, yAux] = rectanch(lenHArm/2,hNode,4);
maze.regx.N5 = [xMain+wNode+lenHArm, yMain+hNode];
maze.regx.N5 = [maze.regx.N5; [NaN NaN]; xAux+wNode+lenHArm, yAux+hNode];
[xAux, yAux] = rectanch(lenHArm/2,hNode,1);
maze.regx.N5 = [maze.regx.N5; [NaN NaN]; xAux+wNode*2+lenHArm, yAux+hNode];

[xMain, yMain] = rectanch(wNode*stickoutMultiplier,hNode*stickoutMultiplier,3);
[xAux, yAux] = rectanch(lenHArm/2,hNode,1);
maze.regx.N6 = [xMain+wNode, yMain+hNode];
maze.regx.N6 = [maze.regx.N6; [NaN NaN]; xAux+wNode, yAux+hNode];

maze.locs.A23 = [xHArm+2*wNode+lenHArm yHArm+h-hNode];
maze.locs.A34 = [xVArm+w-wNode yVArm+hNode];
maze.locs.A45 = [xHArm+2*wNode+lenHArm yHArm];
maze.locs.A56 = [xHArm+wNode yHArm];
maze.locs.A16 = [xVArm yVArm+hNode];
maze.locs.A12 = [xHArm+wNode yHArm+h-hNode];
maze.locs.A25 = [xVArm+wNode+lenHArm yVArm+hNode];

maze.regx.A23 = [xHArm+2*wNode+lenHArm yHArm*stickoutMultiplier+h-hNode];

main = [xVArm*stickoutMultiplier+w-wNode yVArm+hNode];
aux = [[w-wNode; w-wNode-lenHArm/2; w-wNode-lenHArm/2; w-wNode], [2*hNode; 2*hNode; h-2*hNode; h-2*hNode]];
maze.regx.A34 = [main(1:end-1, :); aux; main(end, :)];

maze.regx.A45 = [xHArm+2*wNode+lenHArm yHArm*stickoutMultiplier-2*hNode];

maze.regx.A56 = [w-maze.regx.A45(:, 1), maze.regx.A45(:, 2)];

maze.regx.A16 = [w-maze.regx.A34(:, 1), maze.regx.A34(:, 2)];

maze.regx.A12 = [w-maze.regx.A23(:, 1), maze.regx.A23(:, 2)];

maze.regx.A25 = [xVArm(1)+wNode+lenHArm yVArm(1)+hNode;
    xVArm(2)+wNode+lenHArm yVArm(2)+hNode;
    xVArm(2)+wNode+lenHArm yVArm(2);
    xVArm(2)+wNode+1.5*lenHArm yVArm(2);
    xVArm(3)+wNode+1.5*lenHArm yVArm(3)+2*hNode;
    xVArm(3)+wNode+lenHArm yVArm(3)+2*hNode;
    xVArm(3)+wNode+lenHArm yVArm(3)+hNode;
    xVArm(4)+wNode+lenHArm yVArm(4)+hNode;
    xVArm(4)+wNode+lenHArm yVArm(4)+2*hNode;
    xVArm(4)+wNode+lenHArm/2 yVArm(4)+2*hNode;
    xVArm(4)+wNode+lenHArm/2 yVArm(5);
    xVArm(4)+wNode+lenHArm yVArm(5);
    xVArm(4)+wNode+lenHArm yVArm(5)+hNode];

maze.regx = structfun(@(n) [output(1, @poly2cw, n(:, 1), n(:, 2)), output(2, @poly2cw, n(:, 1), n(:, 2))], maze.regx, 'un', 0);

mazeLocs = fields(maze.locs);
nLocs = length(mazeLocs);
maze.whole.locs = [];
for i = 1:nLocs
   maze.whole.locs = cat(1,maze.whole.locs, [cat(1,maze.locs.(mazeLocs{i})); nan nan]);
end
