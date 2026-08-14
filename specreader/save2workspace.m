function save2workspace(varargin)
% SAVE2WORKSPACE Called by specr and trueref to save current plot data to
%   matlab base workspace.
%
% Copyright 2004, Zhang Jiang

hLine = findall(gca,'Type','line');
if isempty(hLine)
    return;
end

% --- get tags of lines and remove datatipmarkers from hLine
tempHLine = [];;
lineTags = {};
for iLine = 1:length(hLine)
    if ~strcmp(get(hLine(iLine),'Tag'),'DataTipMarker')
        tempHLine(length(tempHLine)+1) = hLine(iLine);
        lineTags{length(lineTags)+1} = get(hLine(iLine),'Tag');
    end
end
hLine = tempHLine';
lineTags = fliplr(lineTags);

for iLine = 1:length(hLine)
    xdata = get(hLine(iLine),'xdata');
    ydata = get(hLine(iLine),'ydata');
    ydataError = getappdata(hLine(iLine),'ydataError');
    tagLine = get(hLine(iLine),'tag');
    assignin('base',tagLine,[xdata' ydata' ydataError]);
end