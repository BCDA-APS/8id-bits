function curvederivative
% CURVEDERIVATIVE Called by specr to calculate curve gradient
%
% Copyright 2009, Zhang Jiang
hFigSpecr = findall(0,'Tag','specr_Fig');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
hLine = findall(gca,'Type','line');

% --- if no curve plotted, return
if isempty(hLine)
    return;
end

colorSpec = varycolor(32);
colorSpec1 = colorSpec(randperm(32),:);
colorSpec2 = colorSpec(randperm(32),:);

% --- calculate gradient
for iLine = length(hLine):-1:1
    x = get(hLine(iLine),'XData');
    y = get(hLine(iLine),'YData');
    lineTag = get(hLine(iLine),'Tag');
    if ~isempty(strmatch('specrscan',lower(lineTag)))
        try
            y1 = gradient(y,x);
            y1 = (y1-min(y1))*(max(y)-min(y))/(max(y1)-min(y1))+min(y);
            hLineDev=line('Parent',hAxes,...
                'XData',x,...
                'YData',y1,...
                'LineStyle',get(hLine(iLine),'LineStyle'),...
                'Color',colorSpec1(iLine,:),...
                'markerfacecolor',colorSpec2(iLine,:),...
                'marker',get(hLine(iLine),'marker'),...
                'markersize',get(hLine(iLine),'markersize'),...
                'Tag',['specrDev',lineTag(10:end)]);
            setappdata(hLineDev,'ydataError',ones(size(y1'))*NaN);            
        catch
        end
    elseif ~isempty(strmatch('specrdev',lower(lineTag)))
        delete(hLine(iLine));        
    end
end

% --- determine legend
curvelegend(hFigSpecr);