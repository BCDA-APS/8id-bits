function truerefimport(varargin)
% TRUEREFIMPORT Called by trueref to import data.
%
% Copyright 2004, Zhang Jiang

hMainFig = findall(0,'Tag','trueref_Fig');
hFigSpecr = findall(0,'Tag','specr_Fig');

settings = getappdata(hFigSpecr,'settings');
openfiles = settings.trueref.openfiles;

% --- figure layout for truerefimport request
figSize = get(hMainFig,'position');
figWidth = 300;
figHeight = 300;
figTruerefimportSize = [figSize(1)+5,(figSize(2)+figSize(4))-figHeight+10,figWidth,figHeight];
hFigTruerefimport = figure(...
    'Units','pixels',...
    'MenuBar','none',...
    'Units','pixels',...
    'Name','Import ASCII Data',...
    'NumberTitle','off',...
    'IntegerHandle','off',...
    'Position',figTruerefimportSize,...
    'HandleVisibility','callback',...
    'Tag','figTruerefimport',...
    'Resize','off',...
    'WindowStyle','modal',...
    'UserData',[]);
panelHeight = 115;
panelSize = [10 figHeight-panelHeight-10 figTruerefimportSize(3)-20 panelHeight];
hPanel = uipanel(...
    'Parent',hFigTruerefimport,...
    'BackgroundColor',get(hFigTruerefimport,'Color'),...
    'Title','Import as',...
    'Units','pixels',...
    'Position',panelSize);
hRadiobutton1 = uicontrol(...
    'Parent',hPanel,...
    'Style','radiobutton',...
    'BackgroundColor',get(hFigTruerefimport,'Color'),...
    'Position',[10 panelSize(4)-40 panelSize(3)-10 25],...
    'String','Reflectivity',...
    'Tag','truerefimport_radiobuttion1',...
    'Value',1,...
    'callback',@radiobuttonFcn);
hRadiobutton2 = uicontrol(...
    'Parent',hPanel,...
    'Style','radiobutton',...
    'BackgroundColor',get(hFigTruerefimport,'Color'),...
    'Position',[10 panelSize(4)-70 panelSize(3)-10 25],...
    'String','Positive longitudinal diffuse (Optional)',...
    'Tag','truerefimport_radiobuttion2',...
    'Value',0,...
    'callback',@radiobuttonFcn);
hRadiobutton3 = uicontrol(...
    'Parent',hPanel,...
    'Style','radiobutton',...
    'BackgroundColor',get(hFigTruerefimport,'Color'),...
    'Position',[10 panelSize(4)-100 panelSize(3)-10 25],...
    'String','Negative longitudinal diffuse (Optional)',...
    'Tag','truerefimport_radiobuttion3',...
    'Value',0,...
    'callback',@radiobuttonFcn);
buttonHeight = 25;
hPushbuttonOK = uicontrol(...
    'Parent',hFigTruerefimport,...
    'Style','pushbutton',...
    'String','Import',...
    'Position',[40 panelSize(2)-40 (figWidth-100)/2 buttonHeight],...
    'callback',{@truerefimport_OKRequestFcn,hMainFig});
hPushbuttonCancel = uicontrol(...
    'Parent',hFigTruerefimport,...
    'Style','pushbutton',...
    'String','Cancel',...
    'Position',[figWidth/2+10 panelSize(2)-40 (figWidth-100)/2 buttonHeight],...
    'Callback',@truerefimport_CloseRequestFcn);
% --- display panel
panelDisplaySize = [10 10 figWidth-20 figHeight-panelHeight-buttonHeight-45];
hPanelDisplay = uipanel(...
    'Parent',hFigTruerefimport,...
    'BackgroundColor',get(hFigTruerefimport,'Color'),...
    'Title','Already imported files',...
    'Units','pixels',...
    'Position',panelDisplaySize);
 hEditDisplayFile = uicontrol(...
    'Parent',hPanelDisplay,...
    'HorizontalAlignment','left',...
    'BackgroundColor','w',...
    'Style','edit',...
    'String',openfiles,...
    'Max',2,...
    'Min',0,...
    'Position',[10 10 panelDisplaySize(3)-15 panelDisplaySize(4)-30],...
    'Enable','on');


%================================================================
% --- radiobutton callback fcn (choose import type)
%================================================================
function radiobuttonFcn(hObject,eventdata)
hFigTruerefimport = gcf;
hRadiobutton1 = findall(hFigTruerefimport,'Tag','truerefimport_radiobuttion1');
hRadiobutton2 = findall(hFigTruerefimport,'Tag','truerefimport_radiobuttion2');
hRadiobutton3 = findall(hFigTruerefimport,'Tag','truerefimport_radiobuttion3');
set(hRadiobutton1,'Value',0);
set(hRadiobutton2,'Value',0);
set(hRadiobutton3,'Value',0);
set(hObject,'Value',1);


%================================================================
% --- truerefimport Ok button callback function
%================================================================
function truerefimport_OKRequestFcn(hObject,eventdata,hMainFig)
% --- get handles of radiobuttons 
hFigTruerefimport = gcf;
hRadiobutton1 = findall(hFigTruerefimport,'Tag','truerefimport_radiobuttion1');
hRadiobutton2 = findall(hFigTruerefimport,'Tag','truerefimport_radiobuttion2');
hRadiobutton3 = findall(hFigTruerefimport,'Tag','truerefimport_radiobuttion3');

hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
openfiles = settings.trueref.openfiles;
openpath  = settings.trueref.openpath;

% --- go the previous save path
prePath = pwd;
if isempty(openpath)
    openpath = pwd;
end
go_openpath_str = ['cd ','''',openpath,''''];
eval(go_openpath_str);
% --- open file dlg
[filename, filepath] = uigetfile( ...
    {'*.dat','Data Files (*.dat)';'*.txt','Text Documents (*.txt)';'*.*','All Files (*.*)'}, ...
    'Select ASCII file','MultiSelect','on');
% If "Cancel" is selected then return
if isequal(filepath,0)
    restorePath(prePath);
    return;
end

settings.trueref.openpath = filepath;
setappdata(hFigSpecr,'settings',settings);
settings = getappdata(hFigSpecr,'settings');
openfiles = settings.trueref.openfiles;
openpath  = settings.trueref.openpath;

% --- get data from main figure or generate new data structure
if isappdata(hMainFig,'data')
    data = getappdata(hMainFig,'data');
else
    data.ref = {};
    data.pos = {};
    data.neg = {};
end
% --- reset merge data and true reflectivity
data.mrg = cell(1,3);
data.tru = [];
setappdata(hMainFig,'data',data);
% --- determine whether multifiles are selected
if iscell(filename)
    nFiles = length(filename);
    for iFiles = 1:nFiles
        file = fullfile(filepath,filename{iFiles});
        try
            imported = load(file);
        catch
            uiwait(msgbox('Invalid ASCII file.','File Open Error','error','modal'));
            rmappdata(hMainFig,'data');
            restorePath(prePath);
            return;
        end
        if size(imported,2)~=3 & size(imported,2)~=2
            uiwait(msgbox('Data should have 2 or 3 columns.','File Open Error','error','modal'));
            rmappdata(hMainFig,'data');
            restorePath(prePath);
            return;
        end
        importFcn(imported,hMainFig,hRadiobutton1,hRadiobutton2,hRadiobutton3,filename{iFiles},filepath);
    end
else
    file = fullfile(filepath,filename);
    try
        imported = load(file);
    catch
        uiwait(msgbox('Invalid ASCII file.','File Open Error','error','modal'));
        restorePath(prePath);
        return;
    end
    if size(imported,2)~=3 & size(imported,2)~=2
        uiwait(msgbox('Data should have 2 or 3 columns.','File Open Error','error','modal'));
        restorePath(prePath);
        return;
    end
    importFcn(imported,hMainFig,hRadiobutton1,hRadiobutton2,hRadiobutton3,filename,filepath);
end
restorePath(prePath);
truerefplot;
%assignin('base','data',data);
delete(gcf);


%================================================================
% --- function called by truerefimport Ok button callback function
%================================================================
function importFcn(imported,hMainFig,hRadiobutton1,hRadiobutton2,hRadiobutton3,filename,filepath);
hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
openfiles = settings.trueref.openfiles;
openpath  = settings.trueref.openpath;
% --- get data
data = getappdata(hMainFig,'data');
% --- determine the third column (error)
if size(imported,2) == 2
    try
        imported(:,3) = sqrt(imported(:,2));
    catch
        imported(:,3) = 0;
    end 
end
switch 1
    case get(hRadiobutton1,'Value')
        data.ref{length(data.ref)+1} = imported;
        openfiles{length(openfiles)+1} = ['Ref:   ',filename];
    case get(hRadiobutton2,'Value')
        data.pos{length(data.pos)+1} = imported;
        openfiles{length(openfiles)+1} = ['Pos:   ',filename];
    case get(hRadiobutton3,'Value')
        data.neg{length(data.neg)+1} = imported;
        openfiles{length(openfiles)+1} = ['Neg:   ',filename];
end
settings.trueref.openfiles = openfiles;
settings.trueref.openpath  = openpath;
setappdata(hMainFig,'data',data);
setappdata(hFigSpecr,'settings',settings);
return;


%================================================================
% --- truerefimport close function
%================================================================
function truerefimport_CloseRequestFcn(hObject,eventdata)
delete(gcf);
return;


%=========================================================
% --- restore to previous path
%=========================================================
function restorePath(prePath)
path_str = ['cd ','''',prePath,''''];
eval(path_str);