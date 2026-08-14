function openspec(varargin)
% OPENSPEC Called by specr to load fids of scans in a selected spec file
%
% Copyright 2004, Zhang Jiang
%
% Line 66: comment out motor names reading sections for bluesky files (2024/10/24) 

hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
% --- get current directory in order to restore after opening spec file
prePath = pwd;
% --- go to the stored directory; 
if ~isempty(file)
    [pathstr,name,ext] = fileparts(file);
    path_str = ['cd ','''',pathstr,''''];
    try 
        eval(path_str);
    catch  % if error, go to current directory
        path_str = ['cd ','''',prePath,''''];
    end
end
% --- open file
[filename, filepath] = uigetfile( ...
    {'*.*','All Files (*.*)'}, ...
    'Select Spec File');
% If "Cancel" is selected then return
if isequal([filename,filepath],[0,0])
    restorePath(prePath);
    return
end
% Otherwise construct the fullfilename and Check and load the file.
file = fullfile(filepath,filename);
[fid,message] = fopen(file,'r');        % open file
if fid == -1                % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'));
    % fclose(fid);
    restorePath(prePath);
    return;
end
scan = '';
scan.motorName = {};
scan.fidPos = [];
scan.head = {};
% set timebar
s=dir(file);
filesize=s.bytes;
posFigSpecr = get(hFigSpecr,'position');
hTimebar = timebar('Please do NOT close while loading scans...','Loading Scans...');
posFigSpecr         = get(hFigSpecr,'position');
posFigTimebarOld    = get(hTimebar,'position');
posFigTimebarNew    = [...
    posFigSpecr(1) + 0.5*(posFigSpecr(3)-posFigTimebarOld(3)),...
    posFigSpecr(2) + 0.5*(posFigSpecr(4)-posFigTimebarOld(4)),...
    posFigTimebarOld(3),...
    posFigTimebarOld(4)];
set(hTimebar,...
    'WindowStyle','modal',...
    'position',posFigTimebarNew);
stopFlag = 0;
while feof(fid) == 0 && stopFlag == 0
    tempfid = ftell(fid);
    scanline = fgetl(fid);

% Commend out for Bluesky files    
%     % --- get motor names
%     if length(scanline) >= 3 && strcmp(scanline(1:2),'#O')
%         % - remove space at the end of the scanline
%         tmpspace = findstr(scanline,' ');
%         while length(scanline) == tmpspace(end)         
%             scanline(end) = '';
%             tmpspace = findstr(scanline,' ');            
%         end
%         % - remove '#Oxx' and space in the front
%         scanline = scanline(tmpspace(1):end);
%         while strcmp(scanline(1),' ')
%             scanline(1) = '';
%         end
%         % - find double space between motor names
%         space = findstr(scanline,'  ');
%         lengthSpace = length(space);
%         for iSpace = lengthSpace:-1:2
%             if space(iSpace) == space(iSpace-1)+1
%                 space(iSpace) = [];
%             end
%         end
%         space = [-1 space length(scanline)+1];
%         % - assign motor names
%         for iCol = 1:(length(space)-1)
%             colHead = scanline(space(iCol)+2:space(iCol+1)-1);
%             while colHead(1) == ' '
%                 colHead(1) = '';
%             end
%             scan.motorName = [scan.motorName,colHead];           
%         end
%     end
    % --- get scan names
    if length(scanline) >= 3 && strcmp(scanline(1:3),'#S ')
        scan.fidPos = [scan.fidPos,tempfid];
        scan.head{length(scan.fidPos)} = scanline;
        try
            timebar(hTimebar,tempfid/filesize);
        catch
            stopFlag = 1;
            break;
        end
    end
end
fclose(fid);     
if stopFlag == 0
    close(hTimebar);
end

if length(scan.head) == 0 | stopFlag == 1
    uiwait(msgbox('Invalid file or no scan in this file','File Open Error','error','modal'));
            restorePath(prePath);   
            resetfigspecr;            
    return;
end
settings.file = file;
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);
settings = getappdata(hFigSpecr,'settings');
restorePath(prePath);
resetfigspecr;
legend off;

%=========================================================
% --- restore to previous path
%=========================================================
function restorePath(prePath)
path_str = ['cd ','''',prePath,''''];
eval(path_str);
