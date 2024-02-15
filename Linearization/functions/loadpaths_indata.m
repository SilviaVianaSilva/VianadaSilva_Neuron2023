function pathData = loadpaths_indata(paths)

inType = class(paths);
switch inType
    case 'char'
        [fn,fp,ext] = fileparts(paths);
        switch ext
            case ''
                pathDirs = {paths};
            case '.txt'
                pathDirs = ReadFileList(paths);
            case '.nvt'
                pathDirs = {[fn '\' fp]};
        end
    case 'cell'
        pathDirs = paths;

    case 'struct'
        pathData = paths;
end


fprintf('Loading %s\n', fullfile(paths,'processedData','indataB.mat'));
load(fullfile(paths,'processedData','indataB.mat'));
disp('Done');
pathData = indata;
%{
nPaths = length(indata);
for p = 1:nPaths
    pathData(p).x = indata(p).x;
    pathData(p).y = indata(p).y;
    %[pathData(p).x, pathData(p).y] = rotatePath(pathData(p).x,pathData(p).y,deg2rad(adj.rotation));
end
%}
%pathData = addVelocityToIndata(pathData);