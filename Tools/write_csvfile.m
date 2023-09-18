function write_csvfile(outfolder, filename, data)
%
%
%
% --------------------------------------------
% Matthias Haberl - UC San Diego 2018

if exist(outfolder, 'dir') ~= 7
mkdir(outfolder);
end

filename_complete = fullfile(outfolder, strcat(filename, '.txt'));
data_written = 0;
dlmwrite(filename_complete,data,'delimiter',';')
data_written = 1;
end