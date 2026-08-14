function invertxaxis(varargin)
hLine = findall(gca,'Type','line');
if isempty(hLine)
    return;
end
for iLine = 1:length(hLine)
    set(hLine(iLine),'XData',-get(hLine(iLine),'XData'));
end

updateparams;