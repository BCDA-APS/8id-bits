function updateparams

hFigSpecr = findall(0,'Tag','specr_Fig');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
hLine = findall(hAxes,'Type','line');
if isempty(hLine)
    return;
end
xdata = get(hLine(end),'XData');
ydata = get(hLine(end),'YData');
scan(:,1) = xdata(:);
scan(:,2) = ydata(:);

%--- calculate Peak, COM, and FWHM and label plot
[scanPeak,scanCOM,scanFWHM] = params(...
    scan);
old_xlabelStr = get(get(hAxes,'XLabel'),'String');
if iscell(old_xlabelStr)
    preStr = old_xlabelStr{1};
else
    preStr = old_xlabelStr;
end
xlabelStr = {preStr;...
    '   ';...
    ['Peak ',num2str(scanPeak.Y),' @ ',num2str(scanPeak.X),...
    ',   COM ',num2str(scanCOM),...
    ',   FWHM ',num2str(scanFWHM.FWHM),' @ ',num2str(scanFWHM.center)]};
set(get(hAxes,'XLabel'),'String',xlabelStr);