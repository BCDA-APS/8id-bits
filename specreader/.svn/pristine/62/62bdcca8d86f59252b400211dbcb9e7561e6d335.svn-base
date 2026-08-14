function truerefmerge(varargin)
% TRUEREFMERGE Called by trueref to load current plots and call mrgker to
%   merge them.
%
% Copyright 2004, Zhang Jiang

hFigSpecr = findall(0,'Tag','specr_Fig');
hFigTrueref = findall(0,'Tag','trueref_Fig');
hAxes = findall(hFigTrueref,'Tag','trueref_Axes');

settings = getappdata(hFigSpecr,'settings');
merge = settings.merge;

% --- return if no data imported
if isappdata(hFigTrueref,'data')
   data = getappdata(hFigTrueref,'data');
else
    return;
end
if isempty(data.ref) & isempty(data.neg) & isempty(data.pos)
    return;
end

% --- merge data by calling truerefmergeker.m
data.mrg = cell(1,3);
try
    if ~isempty(data.ref)
        data.mrg{1} = mrgker(data.ref,merge);
    end
    if ~isempty(data.pos)
        data.mrg{2} = mrgker(data.pos,merge);
    end
    if ~isempty(data.neg)
        data.mrg{3} = mrgker(data.neg,merge);
    end
    setappdata(hFigTrueref,'data',data);
catch
    uiwait(msgbox('Merging failed. Check selected scans.','Merge Error','error','modal'));
    return;
end

% --- determine which curve to plot: all, ref, pos or neg
hLine = findall(hAxes,'Type','line');
if isempty(hLine)
    plot_str = 'All';
else
    % --- remove datatipmarkers from hLine
    tempHLine = [];;
    for iLine = 1:length(hLine)
        if ~strcmp(get(hLine(iLine),'Tag'),'DataTipMarker')
            tempHLine(length(tempHLine)+1) = hLine(iLine);
        end
    end
    hLine = tempHLine';
    tag_str = get(hLine,'Tag');
    if iscell(tag_str)
        tag_str = tag_str';
        for iTag_str = 1:length(tag_str)
            tag_str{iTag_str} = tag_str{iTag_str}(1:3);
        end
        tag_str = cell2mat(tag_str);
    else
        tag_str = tag_str(1:3);
    end
    if ~isempty(findstr(tag_str,'mrg'))
        plot_str = 'All';
    elseif ~isempty(findstr(tag_str,'ref')) &...
            isempty(findstr(tag_str,'pos')) &...
            isempty(findstr(tag_str,'neg'))
        plot_str = 'Ref';
    elseif  isempty(findstr(tag_str,'ref')) &...
           ~isempty(findstr(tag_str,'pos')) &...
            isempty(findstr(tag_str,'neg'))
        plot_str = 'Pos';
    elseif  isempty(findstr(tag_str,'ref')) &...
            isempty(findstr(tag_str,'pos')) &...
           ~isempty(findstr(tag_str,'neg'))
        plot_str = 'Neg';
    else
        plot_str = 'All';
    end
end
switch plot_str
    case 'All'
        iLine_start = 1;
        iLine_end = 3;
    case 'Ref'
        iLine_start = 1;
        iLine_end = 1;
    case 'Pos'
        iLine_start = 2;
        iLine_end = 2;
    case 'Neg'
        iLine_start = 3;
        iLine_end = 3;
end
% --- plot merged data
cla(hAxes);
tag_str = {'mrgReflectivity';'mrgPositiveLongDiffuse';'mrgNegativeLongDiffuse'};
color_str = {'k';'b';'g'};
for iLine = iLine_start:iLine_end
    if ~isempty(data.mrg{iLine})
        hNewLine = line('Parent',hAxes,...
            'XData',data.mrg{iLine}(:,1)',...
            'YData',data.mrg{iLine}(:,2)',...
            'Tag',tag_str{iLine});
        setappdata(hNewLine,'ydataError',data.mrg{iLine}(:,3));
        set(hNewLine,...
            'Color',color_str{iLine},...
            'LineStyle','-',...
            'Marker','o',...
            'MarkerSize',3,...
            'MarkerFaceColor','m');
    end
end
set(hAxes,...
    'XGrid',get(hAxes,'XGrid'),...
    'YGrid',get(hAxes,'YGrid'),...
    'XScale',get(hAxes,'XScale'),...
    'YScale',get(hAxes,'YScale')...
    );

curvelegend(hFigTrueref);
resettoolbar(hFigTrueref);


