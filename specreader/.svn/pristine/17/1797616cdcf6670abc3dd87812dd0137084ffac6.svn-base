function curvesave(varargin)
% CURVESAVE Called by specr and trueref to save current curves to files.
%
% Copyright 2004, Zhang Jiang

hCurrentFig = gcf;
hCurrentAxes = gca;
hLine = findall(gca,'Type','line');

if isempty(hLine)
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
lineTags = fliplr(lineTags);

% --- figure layout for specrsave request
figSize = get(hCurrentFig,'position');
figCenter = [(figSize(1)+figSize(3)/2),...
                  (figSize(1)+figSize(4)/2)];
figSaveSize = [figSize(1)+5,(figSize(2)+figSize(4))-210,240,210];
hFigSave = figure(...
    'Units','pixels',...
    'MenuBar','none',...
    'Units','pixels',...
    'Name','Save to File',...
    'NumberTitle','off',...
    'IntegerHandle','off',...
    'Position',figSaveSize,...
    'HandleVisibility','callback',...
    'Tag','figsave',...
    'Resize','off',...
    'WindowStyle','modal',...
    'UserData',[]);
% --- layeout panel select curve
panelSelectCurveSize = [10 figSaveSize(4)-70 figSaveSize(3)-20 60];
hPanelSelectCurve = uipanel(...
    'Parent',hFigSave,...
    'BackgroundColor',get(hFigSave,'Color'),...
    'Title','Select curve',...
    'Units','pixels',...
    'Position',panelSelectCurveSize);
hPopupmenuSelectCurve = uicontrol(...
    'Parent',hPanelSelectCurve,...
    'Style','popupmenu',...
    'BackgroundColor','w',...
    'Units','pixels',...
    'Position',[10 12.5 panelSelectCurveSize(3)-20 25],...
    'String',lineTags,...
    'HorizontalAlignment','left',...
    'Value',1,...
    'Tag','figsave_PopupmenuSelectCurve');
% --- layout panel saving format
panelSelectFormatSize = [10 figSaveSize(4)-panelSelectCurveSize(4)-90 figSaveSize(3)-20 70];
hPanelSelectFormat = uipanel(...
    'Parent',hFigSave,...
    'BackgroundColor',get(hFigSave,'Color'),...
    'Title','Format',...
    'Units','pixels',...
    'Position',panelSelectFormatSize);
hRadiobutton1 = uicontrol(...
    'Parent',hPanelSelectFormat,...
    'Style','radiobutton',...
    'BackgroundColor',get(hFigSave,'Color'),...
    'Position',[10 panelSelectFormatSize(4)-40 panelSelectFormatSize(3)-20 25],...
    'String','xdata   ydata   error',...
    'Tag','figsave_radiobuttion1',...
    'Value',1,...
    'callback',@radiobutton1Fcn);
hRadiobutton2 = uicontrol(...
    'Parent',hPanelSelectFormat,...
    'Style','radiobutton',...
    'BackgroundColor',get(hFigSave,'Color'),...
    'Position',[10 panelSelectFormatSize(4)-65 panelSelectFormatSize(3)-20 25],...
    'String','xdata   ydata',...
    'Tag','figsave_radiobuttion2',...
    'Value',0,...
    'callback',@radiobutton2Fcn);
% --- layout OK and Cancel buttons
hPushbuttonOK = uicontrol(...
    'Parent',hFigSave,...
    'Style','pushbutton',...
    'String','Next',...
    'Position',[20 20 (figSaveSize(3)-60)/2 25],...
    'callback',{@figsave_OKRequestFcn,hCurrentFig,hCurrentAxes});
hPushbuttonCancel = uicontrol(...
    'Parent',hFigSave,...
    'Style','pushbutton',...
    'String','Cancel',...
    'Position',[figSaveSize(3)/2+10 20 (figSaveSize(3)-60)/2 25],...
    'Callback',@figsave_CloseRequestFcn);


%================================================================
% --- radiobutton1 callback fcn (save without error)
%================================================================
function radiobutton1Fcn(hObject,eventdata)
hFigSave = gcf;
hRadiobutton1 = findall(hFigSave,'Tag','figsave_radiobuttion1');
hRadiobutton2 = findall(hFigSave,'Tag','figsave_radiobuttion2');
set(hRadiobutton1,'Value',1);
set(hRadiobutton2,'Value',0);


%================================================================
% --- radiobutton1 callback fcn (save with error)
%================================================================
function radiobutton2Fcn(hObject,eventdata)
hFigSave = gcf;
hRadiobutton1 = findall(hFigSave,'Tag','figsave_radiobuttion1');
hRadiobutton2 = findall(hFigSave,'Tag','figsave_radiobuttion2');
set(hRadiobutton1,'Value',0);
set(hRadiobutton2,'Value',1);


%================================================================
% --- close function
%================================================================
function figsave_CloseRequestFcn(hObject,eventdata)
delete(gcf);
return;


%================================================================
% --- print to figure
%================================================================
function figsave_OKRequestFcn(hObject,eventdata,hCurrentFig,hCurrentAxes)
% --- determine the curve to save
hFigSave = gcf;
hRadiobutton1 = findall(hFigSave,'Tag','figsave_radiobuttion1');
hRadiobutton2 = findall(hFigSave,'Tag','figsave_radiobuttion2');
hPopupmenuSelectCurve = findobj(hFigSave,'Tag','figsave_PopupmenuSelectCurve');
curveStrList = get(hPopupmenuSelectCurve,'String');
curveStr = curveStrList{get(hPopupmenuSelectCurve,'Value')};
hLine = findall(hCurrentAxes,'Type','line');
hLine2Save = findall(hLine,'Tag',curveStr);

% --- get setting parameters
hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
savepath = settings.savepath;

% --- get data from curve
xdata = get(hLine2Save,'xdata');
ydata = get(hLine2Save,'ydata');
ydataError = getappdata(hLine2Save,'ydataError');
if get(hRadiobutton1,'Value') == 1      % determine saving w/t error
    data2save = [xdata' ydata' ydataError];
else
    data2save = [xdata' ydata'];
end
% --- go the previous save path
prePath = pwd;
if isempty(savepath)
    savepath = pwd;
end

try
    go_savepath_str = ['cd ','''',savepath,''''];
    eval(go_savepath_str);
catch
    go_savepath_str = ['cd ','''',pwd,''''];
    eval(go_savepath_str);
end
% --- user choose a file to save
filter={'*.dat','Data File (*.dat)';...
        '*.txt','Text Documents (*.txt)';...
        '*.*','All Files (*.*)'};
[filename, pathname,filterIndex] = uiputfile(filter,'Save Data As',[curveStr,'.dat']);	
% If 'Cancel' was selected then return
if isequal([filename,pathname],[0,0])
    restorePath(prePath);
    return
else
    % Construct the full path and save
    File = fullfile(pathname,filename);
    switch filterIndex
        case 1
            if ~ strcmp(File(end-3:end),'.dat')
                File=[File '.dat'];
            end
        case 2
            if ~ strcmp(File(end-3:end),'.txt')
                File=[File '.txt'];
            end
        case 3
    end
    go_savepath_str = ['cd ','''',pathname,''''];    % change path
    eval(go_savepath_str);
    save_str=['save ',filename,' data2save -ascii -tabs'];
    eval(save_str);
end
restorePath(prePath);
settings.savepath = pathname;
setappdata(hFigSpecr,'settings',settings);
settings = getappdata(hFigSpecr,'settings');
delete(hFigSave);


%=========================================================
% --- restore to previous path
%=========================================================
function restorePath(prePath)
path_str = ['cd ','''',prePath,''''];
eval(path_str);
