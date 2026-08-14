function plotstyle(varargin)
hPopupmenuPlotStyle = gcbo;
hAxes = gca;
if isempty(hAxes) | isempty(hPopupmenuPlotStyle)
    return;
end

switch get(hPopupmenuPlotStyle,'Value')
    case 1          % linear
        set(hAxes,'XScale','linear','YScale','linear');
    case 2          % logx
        set(hAxes,'XScale','log','YScale','linear');
    case 3          % logy
        set(hAxes,'XScale','linear','YScale','log');
    case 4          % logxy
        set(hAxes,'XScale','log','YScale','log');
end
