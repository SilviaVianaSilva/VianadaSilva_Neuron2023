% Import all data from the NeuraLynx video tracker data file
function handles = readVideoData(file)

fieldSelect = [1,1,1,1,1,1];
getHeader = 0;
extractMode = 1;

% Get the data
[handles.post,handles.posx,handles.posy,handles.angles,handles.targets,handles.points] =...
    Nlx2MatVT(file,fieldSelect,getHeader,extractMode);

% Convert data
handles.targets = uint32(handles.targets);