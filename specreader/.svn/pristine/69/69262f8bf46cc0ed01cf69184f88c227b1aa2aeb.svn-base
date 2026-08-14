function scanmerge(varargin)
% SCANMERGE Called by specr to read scans from current plot to call mrgker
%   to merge scans.
%
% Copyright 2004, Zhang Jiang

hFigSpecr = findall(0,'Tag','specr_Fig');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
hLine = findall(hAxes,'Type','line');

% --- get merge mode
settings = getappdata(hFigSpecr,'settings');
merge = settings.merge;

% --- if no or only one curve plotted, return
if isempty(hLine) | length(hLine) <= 1
    return;
end

% --- if there are derivative curves, return;
for iLine = 1:length(hLine)
    if ~isempty(strmatch('specrdev',lower(get(hLine(iLine),'Tag'))))
        uiwait(msgbox('Merging failed. Check selected scans.','Merge Error','error','modal'));        
        return; 
    end;
end

% --- get tags of lines and remove datatipmarkers from hLine
tempHLine = [];
lineTags = {};
for iLine = 1:length(hLine)
    if ~strcmp(get(hLine(iLine),'Tag'),'DataTipMarker')
        tempHLine(length(tempHLine)+1) = hLine(iLine);
        lineTags{length(lineTags)+1} = get(hLine(iLine),'Tag');
    end
end
hLine = tempHLine';
lineTags = fliplr(lineTags);

% --- extract x,y data, yerror, scan tags, xlable, ylabel, title from each plot
scanTag = {};
scanData = {};
scanIndex = [];
for iLine = 1:length(hLine)
    xdata = get(hLine(iLine),'xdata');
    ydata = get(hLine(iLine),'ydata');
    ydataError = getappdata(hLine(iLine),'ydataError');
    scanTag{iLine} = get(hLine(iLine),'tag');
    scanIndex = [scanIndex str2num(scanTag{iLine}(10:end))];
    scanData{iLine} = [xdata' ydata' ydataError];
end
xlabel_str = get(get(hAxes,'XLabel'),'String');
ylabel_str = get(get(hAxes,'YLabel'),'String');
scanIndex = sort(scanIndex);
% construct scan str
b={};
b{1,1} = scanIndex(1);
for ii=2:length(scanIndex)
    if scanIndex(ii)-scanIndex(ii-1) == 1
        b{end,1} = [b{end,1},scanIndex(ii)];
    else
        b{end+1,1} = scanIndex(ii);
    end
end
scanStr = '';
for ii=1:length(b)
    if length(b{ii})>1
       tmp_str = [num2str(b{ii}(1)),'-',num2str(b{ii}(end))];
    else
       tmp_str = num2str(b{ii});
    end
    scanStr = [scanStr,' ',tmp_str];    
end
scanStr(1) = '';
% generate title str
[~,filename,~] = fileparts(settings.file);
title_str = {...
    ['File: ',titlestr(filename),',   Merged scans ',scanStr];...
    ['']};
%title_str = ['Merged scans ',num2str(scanIndex)];

% --- if merging failed do nothing
try
    mergeData = mrgker(scanData,merge);  % call mrgker function to merge scanData
catch
    uiwait(msgbox('Merging failed. Check selected scans.','Merge Error','error','modal'));
    return;
end

% --- clear previous plot and plot merged data
cla(hAxes);
if length(scanIndex) == 1
    mergeLineTag    = ['specrMerge',num2str(scanIndex)];
else
    mergeLineTag    = ['specrMerge',num2str(min(scanIndex)),'To',num2str(max(scanIndex))];
end
hNewLine = line('Parent',hAxes,...
        'XData',mergeData(:,1)',...
        'YData',mergeData(:,2)',...
        'Tag',mergeLineTag);
setappdata(hNewLine,'ydataError',mergeData(:,3));
set(hNewLine,...
    'Color','b',...
    'LineStyle','-',...
    'Marker','o',...
    'MarkerSize',5,...
    'MarkerFaceColor','m');

% --- xlabel,ylabel and title
if iscell(xlabel_str)
    xlabel_str = xlabel_str{1};
end
set(hAxes,...
    'XGrid',get(hAxes,'XGrid'),...
    'YGrid',get(hAxes,'YGrid'),...
    'XScale',get(hAxes,'XScale'),...
    'YScale',get(hAxes,'YScale'),...
    'XLabel',text('String',xlabel_str),...
    'YLabel',text('String',ylabel_str),...
    'Title',text('String',title_str)...
    );
curvelegend(hFigSpecr);
resettoolbar(hFigSpecr);

updateparams;