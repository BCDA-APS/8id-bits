function refreshscanplot(varargin)
% REFRESHSCANPLOT Called by specr to reload scan(s) and plot.
%
% Copyright 2004, Zhang Jiang

hFigSpecr = gcf;
hPopupmenuX = findall(hFigSpecr,'Tag','specr_PopupmenuX');
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
hPopupmenuPlotStyle = findall(hFigSpecr,'Tag','specr_PopupmenuPlotStyle');
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
if isempty(file) | isempty(scan) | ~isfield(scan,'selectionIndex') | isempty(scan.selectionIndex)
    return;
end

% check dynamically if there are more scans added to spec file
[fid,message] = fopen(file);
fseek(fid,scan.fidPos(end),'bof');  % move fid to end of the last stored scan
scanline = fgetl(fid);      % skip the currentline (head of the last stored scan);
while feof(fid) == 0
    tempfid = ftell(fid);
    scanline = fgetl(fid);
    if length(scanline) >= 3 && strcmp(scanline(1:3),'#S ')
        scan.fidPos = [scan.fidPos,tempfid];
        scan.head{length(scan.fidPos)} = scanline;
    end
end
fclose(fid);
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);
settings = getappdata(hFigSpecr,'settings');
selection = scan.selectionIndex;
scan.selection = cell(length(selection),1);
[fid,message] = fopen(file);
for iSelection = 1:length(selection)
    fseek(fid,scan.fidPos(selection(iSelection)),'bof');
    % --- read time of scan #D
    scanline = fgetl(fid);
    while ~strcmp(scanline(1:2),'#D')
        scanline = fgetl(fid);
    end
    scan.selection{iSelection}.time = scanline(4:end);  
    % --- read current motor positions
    scan.selection{iSelection}.motorPos = [];
    while ~strcmp(scanline(1:2),'#P')
        scanline = fgetl(fid);
    end
    while length(scan.selection{iSelection}.motorPos) < length(scan.motorName)
        tmpspace = findstr(scanline,' ');
        scanline = scanline(tmpspace(1):end);
        scan.selection{iSelection}.motorPos = [scan.selection{iSelection}.motorPos, str2num(scanline)];
        scanline = fgetl(fid);
    end    
    % --- read number of scan columns #N ...
    while ~strcmp(scanline(1:2),'#N')   
        scanline = fgetl(fid);
    end
    scan.selection{iSelection}.length = str2num(scanline(3:end));
    % --- read head of scan columns #L ...
    while ~strcmp(scanline(1:2),'#L')
        scanline = fgetl(fid);
    end
    scanline = scanline(4:end);
    space = findstr(scanline,'  ');
    lengthSpace = length(space);
    for iSpace = lengthSpace:-1:2
        if space(iSpace) == space(iSpace-1)+1
            space(iSpace) = [];
        end
    end
    space = [-1 space length(scanline)+1];
    scan.selection{iSelection}.colHead = ...
        cell(1,scan.selection{iSelection}.length);
    for iCol = 1:scan.selection{iSelection}.length
        colHead = scanline(space(iCol)+2:space(iCol+1)-1);
        while colHead(1) == ' '
            colHead(1) = '';
        end
        scan.selection{iSelection}.colHead{iCol} = colHead;
    end
    
    % --- read scan data
    colData = [];
    str_scanline = num2str(scanline);
    while ~strcmp(str_scanline,'') & ~strcmp(str_scanline,'-1')
        fidPos = ftell(fid);
        scanline = fgetl(fid);
        str_scanline = num2str(scanline);
        while ~strcmp(str_scanline,'') & ~strcmp(str_scanline,'-1')...
                & ~strcmp(str_scanline(1:1),'#')
            fseek(fid,fidPos,'bof');
            colData = [colData fscanf(fid,'%g',[scan.selection{iSelection}.length,1])];
            scanline = fgetl(fid);
            fidPos = ftell(fid);
            scanline = fgetl(fid);
            str_scanline = num2str(scanline);
        end
    end
    scan.selection{iSelection}.colData = colData';
end
fclose(fid);
% --- if not the same type scan selected, error and return;
if length(selection) >= 2
    for iSelection = 2:length(selection)
        if scan.selection{iSelection}.length == scan.selection{iSelection-1}.length
            colHeadCmp = strcmp(scan.selection{iSelection}.colHead,...
                scan.selection{iSelection-1}.colHead);
        else
            colHeadCmp = 0;
        end
        if ~isempty(find(colHeadCmp == 0))
            uiwait(msgbox('Multi-selections have to be scans of the same type.',...
                'Select Scan Error','error','modal'));
            return;
        end
    end
end
scan.selectionIndex = selection;        % store the selected scan index
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);
settings = getappdata(hFigSpecr,'settings');
resettoolbar(hFigSpecr);
try
    scanplot;
catch
end
zoom out;