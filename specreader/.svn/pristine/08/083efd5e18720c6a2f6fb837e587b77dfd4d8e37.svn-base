function selectscan(varargin)
% SELECTSCAN Called by specr to open scan listbox and load scans to plot.
%
% Copyright 2004, Zhang Jiang

hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
if isempty(file) | isempty(scan.fidPos)
    return;
end

% --- update spec file
updatespecscan(hFigSpecr);

% --- read scan data
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
[filepath,filename,fileext] = fileparts(file);
try             % highlight the current selection or last scan
    highlightScanIndex = sort(scan.selectionIndex);
catch
    highlightScanIndex = length(scan.fidPos);
end
[selection,ok] = listdlg(...
    'PromptString',['File: ',filename,fileext],...
    'SelectionMode','multiple',...
    'Name','Select Scan',...
    'OKString','Plot',...
    'ListSize',[300 400],...
    'ListString',scan.head,...
    'InitialValue',highlightScanIndex);
if isempty(selection) | ok == 0
    return;
end
% read scan
scan = readspecscan(file,scan,selection);

% --- post selection scan
postselectscan(hFigSpecr,scan,selection);