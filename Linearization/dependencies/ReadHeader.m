function H = ReadHeader(fp)
% H = ReadHeader(fp)
%  Reads NSMA header, leaves file-read-location at end of header
%  INPUT:

%      fid -- file-pointer (i.e. not filename)
%  OUTPUT:
%      H -- cell array.  Each entry is one line from the NSMA header
% Now works for files with no header.
% ADR 1997
% version L4.1
% status: PROMOTED
% v4.1 17 nov 98 now works for files sans header
%
% May 2010, modified by Emily to work with new Matlab
%---------------

% Get keys
beginheader = '%%BEGINHEADER';
endheader = '%%ENDHEADER';

iH = 1; H = {};
curfpos = ftell(fp);

% look for beginheader
headerLine = my_fgetl(fp);
if strcmp(headerLine, beginheader)
    H{1} = headerLine;
    while ~feof(fp) & ~strcmp(headerLine, endheader)
        headerLine = my_fgetl(fp);
        iH = iH+1;
        H{iH} = headerLine;
    end
else % no header
    fseek(fp, curfpos, 'bof');
end

function tline = my_fgetl(fid)

try
    [tline,lt] = fgets(fid);
    tline = tline(1:end-length(lt));
    fseek(fid, -(length(lt)-1), 'cof');
    
    if isempty(tline)
        tline = '';
    end
    
catch exception
    if nargin ~= 1
        error (nargchk(1,1,nargin,'struct'))
    end
    throw(exception);
end


% function H = ReadHeader(fp)
% 
% 
% % H = ReadHeader(fp)
% %
% %  Reads NSMA header, leaves file-read-location at end of header
% %
% %  INPUT: 
% %
% %      fid -- file-pointer (i.e. not filename)
% %  OUTPUT: 
% %
% %      H -- cell array.  Each entry is one line from the NSMA header
% %
% %
% %
% % Now works for files with no header.
% %
% % ADR 1997
% %
% % version L4.1
% %
% % status: PROMOTED
% %
% %
% %
% % v4.1 17 nov 98 now works for files sans header
% %
% 
% % ******************* THIS CODE WAS TAKEN FROM MCLUST. DON'T MODIFY IT OR
% % THERE WILL BE INTERFERING VERSIONS. IF YOU WANT SOMETHING THAT DOES
% % SOMETHING SIMILAR, COPY THIS AND CALL IT SOMETHING ELSE! **************
% 
% %---------------
% 
% % Get keys
% 
% beginheader = '%%BEGINHEADER';
% endheader = '%%ENDHEADER';
% 
% iH = 1; H = {};
% 
% curfpos = ftell(fp);
% 
% 
% %--------------
% 
% % go
% 
% 
% 
% % look for beginheader
% 
% headerLine = fgetl(fp);
% 
% if strcmp(headerLine, beginheader)
% 
%    H{1} = headerLine;
% 
%    while ~feof(fp) & ~strcmp(headerLine, endheader)     
%       headerLine = fgetl(fp);
%       iH = iH+1;
% 
%       H{iH} = headerLine;
% 
%    end
% 
% else % no header
% 
%    fseek(fp, curfpos, 'bof');
% end
