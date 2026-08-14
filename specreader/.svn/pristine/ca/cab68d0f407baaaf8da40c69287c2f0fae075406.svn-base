function curvenormalize(varargin)
% CURVENORMALIZE Called by specr and trueref
%
% Copyright 2004, Zhang Jiang

hLine = findall(gca,'Type','line');
% --- if no curve plotted, return
if isempty(hLine)
    return;
end

% --- dlg to input normalization factor
dlgAnswer=inputdlg({'Enter normalization factor: '},'Normalization',1);
if isempty(dlgAnswer)
    return;
end
normFactor = str2double(dlgAnswer{1});
if isnan(normFactor) | isempty(normFactor) | normFactor <=0
    uiwait(msgbox('Invalid argument!','Normalization Error','error','modal'));
    return;
end
for iLine = 1:length(hLine)
    set(hLine(iLine),'YData',get(hLine(iLine),'YData')/normFactor);
    ydataError =  getappdata(hLine(iLine),'ydataError');
    ydataError = ydataError/normFactor;
    setappdata(hLine(iLine),'ydataError',ydataError);
end


