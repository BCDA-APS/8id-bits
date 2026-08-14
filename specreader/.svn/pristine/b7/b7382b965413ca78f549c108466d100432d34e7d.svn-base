function scanplot(varargin)
% SCANPLOT Called by specr to plot scans.
%
% Copyright 2009, Zhang Jiang
%
% Line 95: add auto axis for bluesky files due to nan xlims. The changes are good for regular spec files too (2024/10/24). 

hFigSpecr = findall(0,'Tag','specr_Fig');
hPopupmenuX = findall(hFigSpecr,'Tag','specr_PopupmenuX');
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
hPopupmenuPlotStyle = findall(hFigSpecr,'Tag','specr_PopupmenuPlotStyle');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
settings = getappdata(hFigSpecr,'settings');

file = settings.file;
scan = settings.scan;

if isempty(file) | isempty(scan) | ~isfield(scan,'selection') | isempty(scan.selection{1})
    return;
end

% --- set x, y axis labels
popupmenuXString = get(hPopupmenuX,'String');
popupmenuYString = get(hPopupmenuY,'String');
popupmenuXValue = get(hPopupmenuX,'Value');
popupmenuYValue = get(hPopupmenuY,'Value');
set(get(hAxes,'XLabel'),'String',popupmenuXString(popupmenuXValue));
set(get(hAxes,'YLabel'),'String',popupmenuYString(popupmenuYValue));


cla(hAxes);
hLine = [];
colorSpec = {'b';'r';'g';'c';'m';'k'};
markerSpec = {'o';'s';'d';'^';'v';'<';'>'};
% markerFaceColorSpec = varycolor(32);
% markerFaceColorSpec = markerFaceColorSpec(randperm(32),:);
markerFaceColorSpec = {'m';'b';'k';'r';'c';'g'};
for iSelection = 1:length(scan.selection)
    % --- if ydata has negative data,then ydataError = 0
    if isempty(find( scan.selection{iSelection}.colData(:,popupmenuYValue)<0 ))
        % --- statistical error
        ydataError = ...
            sqrt(scan.selection{iSelection}.colData(:,popupmenuYValue));
    else
        ydataError = ...
            zeros(length(scan.selection{iSelection}.colData(:,popupmenuYValue)),1);
    end

    hLine(iSelection) = line('Parent',hAxes,...
        'XData',scan.selection{iSelection}.colData(:,popupmenuXValue),...
        'YData',scan.selection{iSelection}.colData(:,popupmenuYValue),...
        'Tag',['specrScan',num2str(scan.selectionNumber(iSelection))]);
    set(hLine(iSelection),...
        'color',colorSpec{mod(iSelection-1,length(colorSpec))+1},...
        'marker',markerSpec{mod(iSelection-1,length(markerSpec))+1},...
        'markerfacecolor',markerFaceColorSpec{mod(iSelection-1,length(markerFaceColorSpec))+1});
%     set(hLine(iSelection),...
%         'color',colorSpec{mod(iSelection-1,length(colorSpec))+1},...
%         'marker',markerSpec{mod(iSelection-1,length(markerSpec))+1},...
%         'markerfacecolor',markerFaceColorSpec(mod(iSelection-1,length(markerFaceColorSpec))+1,:));
     setappdata(hLine(iSelection),'ydataError',ydataError);
end
set(hLine,...
    'LineStyle','-',...
    'MarkerSize',5);

% set(hLine,...
%     'Color','b',...
%     'LineStyle','-',...
%     'Marker','o',...
%     'MarkerSize',5,...
%     'MarkerFaceColor','m');

% --- determine whether to fix xlim for a monitor scan
% fix xlim only when the first column is selected and it is not a time scan
% nor a hklscan
if ~isempty(timerfindall('Tag','timerMonitor')) ...
    && settings.monitorErasemode == 1 ...
    && popupmenuXValue == 1 ...
    && isempty(findstr(lower(scan.head{end}),' timescan'))  ...
    && isempty(findstr(lower(scan.head{end}),' xpcsscan')) ...
    && isempty(findstr(lower(scan.head{end}),' hklscan')) 
    a = scan.head{end};
    b = findstr(a,' ');
    c = diff(b);
    for ii = length(c):-1:1
        if c(ii)==1
            a(b(ii)) = [];
        end
    end
    b = findstr(a,' ');
    xlims = sort([str2double(a(b(4)+1:b(5)-1)), str2double(a(b(5)+1:b(6)-1))]);
    xlims = [xlims(1)-(xlims(2)-xlims(1))*0.05-eps, xlims(2)+(xlims(2)-xlims(1))*0.05+eps];

    % skip xlims for bluesky files due to NaN for xlims (2024/10/24)
    if any(isnan(xlims))
        set(hAxes,'xlimmode','auto');
    else
       set(hAxes,'xlim',xlims); 
    end

else
    set(hAxes,'xlimmode','auto');
end

% --- set title
[filepath,filename,fileext] = fileparts(file);
% construct scan str
b={};
b{1,1} = scan.selectionNumber(1);
for ii=2:length(scan.selectionNumber)
    if scan.selectionNumber(ii)-scan.selectionNumber(ii-1) == 1
        b{end,1} = [b{end,1},scan.selectionNumber(ii)];
    else
        b{end+1,1} = scan.selectionNumber(ii);
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
% construct title str
title_str = {...
    ['File: ',titlestr(filename),',   Scan ',scanStr,',   ',scan.selection{end}.time];...
    [scan.head{scan.selectionIndex(end)}(length(num2str(scan.selectionNumber(end)))+5:end)]...
    };
set(get(hAxes,'Title'),'String',title_str);

%--- calculate Peak, COM, and FWHM and label plot
[scanPeak,scanCOM,scanFWHM] = params(...
    scan.selection{end}.colData(:,[popupmenuXValue,popupmenuYValue]));
xlabelStr = {scan.selection{end}.colHead{get(hPopupmenuX,'Value')};...
    '   ';...
    ['Peak ',num2str(scanPeak.Y),' @ ',num2str(scanPeak.X),...
    ',   COM ',num2str(scanCOM),...
    ',   FWHM ',num2str(scanFWHM.FWHM),' @ ',num2str(scanFWHM.center)]};
set(get(hAxes,'XLabel'),'String',xlabelStr);

% --- determine legend
curvelegend(hFigSpecr);