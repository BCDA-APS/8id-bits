function invertyaxis(varargin)
% Copyright 2012, Zhang Jiang
hLine = findall(gca,'Type','line');
if isempty(hLine)
    return;
end
for iLine = 1:length(hLine)
        ydata = get(hLine(iLine),'YData');
        y_max = max(ydata);
        y_min = min(ydata);
        set(hLine(iLine),'YData',y_max+y_min-ydata);
end

updateparams;


