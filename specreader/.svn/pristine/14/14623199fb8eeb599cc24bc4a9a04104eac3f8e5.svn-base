function subtractbkg(varargin)
% SUBTRACTBKG GUI to subtract constant background; called by specr.
% 
% Copyright 2006, Zhang Jiang


hFigSpecr = findall(0,'Tag','specr_Fig');
hAxes = findall(hFigSpecr,'Tag','specr_Axes');
hLine = findall(hAxes,'Type','line');
if isempty(hLine)
    return;
end
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
nOfScan = length(scan.selectionIndex);
if isempty(file) | isempty(scan) | ~isfield(scan,'selection') | isempty(scan.selection{1})
    return;
end


% --- get tags of lines and remove datatipmarkers from hLine
tempHLine = [];;
lineTags = {};
for iLine = 1:length(hLine)
    if ~strcmp(get(hLine(iLine),'Tag'),'DataTipMarker')
        tempHLine(length(tempHLine)+1) = hLine(iLine);
        lineTags{length(lineTags)+1} = get(hLine(iLine),'Tag');
    end
end
hLine = tempHLine';
hLine = flipud(hLine);
lineTags = fliplr(lineTags);


% --- check whether the current curve is not specScan. If not, return;
for iLine = 1:length(hLine)
    if ~strcmp(lineTags{iLine}(1:9),'specrScan');
        uiwait(msgbox('Background subtraction does not apply.','Background Subtraction Error','error','modal'));
        return;
    end
end


% --- dlg
prompt  = cell(nOfScan,1);
initial = cell(nOfScan,1);
for iOfScan = 1:nOfScan
    prompt{iOfScan}     = scan.head{scan.selectionIndex(iOfScan)};
    initial{iOfScan}    = num2str(0);
end
dlgAnswer = inputdlg(prompt,'Constant Background Subtraction',1,initial,'on');
if isempty(dlgAnswer)
    return;
end
bkg = zeros(1,length(nOfScan));
for iOfScan = 1:nOfScan
    temp = str2double(dlgAnswer{iOfScan});
    if isnan(temp)
        uiwait(msgbox(['Invalid background for scan ',num2str(scan.selectionNumber(iOfScan))],'Background Subtraction Error','error','modal'));        
        return;
    end
    bkg(iOfScan) = str2double(dlgAnswer{iOfScan});
end


% --- subtract background and update the plot
for iLine = 1:length(hLine)
    xdata = get(hLine(iLine),'XDATA');
    ydata = get(hLine(iLine),'YDATA');
    ydataRelativeError = getappdata(hLine(iLine),'ydataError')./ydata';
    ydata = ydata - bkg(iLine);
    ydataError  = ydataRelativeError.*ydata';
    negIndex    = find(ydata<=0);
    xdata(negIndex)         = [];
    ydata(negIndex)         = [];
    ydataError(negIndex)    = [];
    set(hLine(iLine),'XDATA',xdata);
    set(hLine(iLine),'YDATA',ydata);
    setappdata(hLine(iLine),'ydataError',ydataError);
end

updateparams;
