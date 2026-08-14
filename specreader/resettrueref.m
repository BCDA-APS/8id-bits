function resettrueref(varargin)
hFigTrueref = findall(0,'Tag','trueref_Fig');
hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
if isempty(hFigTrueref)
    return;
end
if isappdata(hFigTrueref,'data');
    rmappdata(hFigTrueref,'data');  % delete figure's children 'data'
    settings.trueref.openfiles = {};
    setappdata(hFigSpecr,'settings',settings);
end
hAxes = findall(hFigTrueref,'Tag','trueref_Axes');
if ~isempty(hAxes)
    cla(hAxes);
end
legend off;
resettoolbar(hFigTrueref);