function curvesmooth(varargin)
% CURVESMOOTH Called by specr to smooth current scans
%
% Copyright 2009, Zhang Jiang

hLine = findall(gca,'Type','line');
% --- if no curve plotted, return
if isempty(hLine)
    return;
end
% --- dlg to input normalization factor
for iLine = 1:length(hLine)
    x = get(hLine(iLine),'XData');
    y = get(hLine(iLine),'YData');
    try
        y1 = smooth(x,y);
        set(hLine(iLine),'YData',y1);
    catch
    end
end

updateparams;
