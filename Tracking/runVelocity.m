
for i = [1450] 

fprintf('\nRunning animal i %d',i);     
    
close all                                   
clearvars -except i
load('C:\Users\Silvia\Dropbox\SilviaProjectCode\+Figure8DataOrganization\figure8Data.mat')
Target_size_maze = [52 72]; % Insert here the total length of the maze
[I, scaleB] = writeVelocity.getPreprocessedIndata(path{i},bSess{i},diode(i),Target_size_maze,1); 
scaleB

% and then we run through sleep
Target_size_box = [31.5 23]; % Insert here the total length of the box, will not be used for scaling though
[I, scaleS]= writeVelocity.getPreprocessedIndata(path{i},sSess{i},diode(i),Target_size_box,0, 'scaleB',scaleB);
%close all
clear I;

end
