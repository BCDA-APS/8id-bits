function curvelegend(hCurrentFig)
% CURVELENGEND Called by specr and trueref
%
% Copyright 2004, Zhang Jiang

hToolbarLegend = findall(hCurrentFig,'Tag','toolbarLegend');
hAxes = findall(hCurrentFig,'Type','axes');
hAxes = hAxes(end);
hLine = findall(hAxes,'Type','line');

if isempty(hLine)
    set(hToolbarLegend,'state','off');
    legend off;
    return;
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

if strcmp(get(hToolbarLegend,'state'),'off')
    set(hToolbarLegend,'state','off');
    delete(findall(hCurrentFig,'Tag','legend'));
    return;
end
legend_str = {};

for iLine = 1:length(hLine)
    iLineTag = get(hLine(iLine),'Tag');
    legend_str{iLine} = iLineTag;
end
legend(hAxes,fliplr(legend_str));

