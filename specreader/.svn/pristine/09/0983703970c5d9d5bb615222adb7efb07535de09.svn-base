function resetfigspecr(varargin)
% Reset all the plots in specr
hFigSpecr = findall(0,'Tag','specr_Fig');
hPopupmenuX = findall(hFigSpecr,'Tag','specr_PopupmenuX');
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
hPopupmenuPlotStyle = findall(hFigSpecr,'Tag','specr_PopupmenuPlotStyle');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
settings = getappdata(hFigSpecr,'settings');

resettoolbar(hFigSpecr);
cla(hAxes);
set(hPopupmenuX,'String','X Axis','Value',1);
set(hPopupmenuY,'String','Y Axis','Value',1);

file = settings.file;
[filepath,filename,fileext] = fileparts(file);
set(hAxes,'Title',text('String',['File: ',titlestr([filename,fileext])]));
set(hAxes,'XLabel',text('String',''));
set(hAxes,'YLabel',text('String',''));
set(hFigSpecr,'Name',['Spec Reader - ',filename,fileext]);

hTimerMonitor   = timerfindall('Tag','timerMonitor');
if ~isempty(hTimerMonitor)
    if isvalid(hTimerMonitor)
        stop(hTimerMonitor);
    end
    delete(hTimerMonitor);
end