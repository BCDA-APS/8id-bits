function convert2theta2qz(varargin)
hFigSpecr = findall(0,'Tag','specr_Fig');
hPopupmenuX = findall(hFigSpecr,'Tag','specr_PopupmenuX');
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
hPopupmenuPlotStyle = findall(hFigSpecr,'Tag','specr_PopupmenuPlotStyle');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
hLine = findall(hAxes,'Type','line');
settings = getappdata(hFigSpecr,'settings');
wavelength = settings.wavelength;
file = settings.file;
scan = settings.scan;
if isempty(file) | isempty(scan) | ~isfield(scan,'selectionIndex') | isempty(scan.selectionIndex)
    return;
end
for iLine = 1:length(hLine)
    set(hLine(iLine),'XDATA',4*pi/wavelength*sin(1/2*pi/180*abs(get(hLine(iLine),'XDATA'))));
end
set(hAxes,'XLabel',text('String','Qz(A^{-1})'));
