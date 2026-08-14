function monitorscan
% OPENSPEC Called by specr to load fids of scans in a selected spec file
%
% Copyright 2005, Zhang Jiang

hFigSpecr = findall(0,'Tag','specr_Fig');
hToolbarMonitor = findall(hFigSpecr,'Tag','toolbarMonitor');
hTimerMonitor   = timerfindall('Tag','timerMonitor');

% --- get setttings
settings = getappdata(hFigSpecr,'settings');
settings.monitorAutoPeriod = settings.monitorPeriod;
setappdata(hFigSpecr,'settings',settings);
file = settings.file;
scan = settings.scan;

% --- get initial 'enable' property of menus,toolbars and pushbuttons
hEnableGroup  = [...
    findall(hFigSpecr,'Tag','specr_MenuTools');...
    findall(hFigSpecr,'Tag','specr_MenuFile');...
    findall(hFigSpecr,'Tag','toolbarOpen');...
    findall(hFigSpecr,'Tag','toolbarSave');...
    findall(hFigSpecr,'Tag','toolbarPrint');...
    findall(hFigSpecr,'Tag','toolbarEditPlot');...
    findall(hFigSpecr,'Tag','toolbarZoom');...
    findall(hFigSpecr,'Tag','toolbarDataCursor');...
    findall(hFigSpecr,'Tag','toolbarLegend');...
    findall(hFigSpecr,'Tag','toolbarPlottoolsOff');...
    findall(hFigSpecr,'Tag','toolbarPlottoolsOn');...
    findall(hFigSpecr,'Tag','toolbarScanPre');...    
    findall(hFigSpecr,'Tag','toolbarScanNext');...    
    findall(hFigSpecr,'Tag','toolbarInvert');...
    findall(hFigSpecr,'Tag','toolbarInvertY');...
    findall(hFigSpecr,'Tag','toolbarSubtractbkg');...
    findall(hFigSpecr,'Tag','toolbarSmooth');...
    findall(hFigSpecr,'Tag','toolbarDev');...
    findall(hFigSpecr,'Tag','toolbarShowScan');...
    findall(hFigSpecr,'Tag','toolbarShowMotor');...
    findall(hFigSpecr,'Tag','toolbarMerge');...
    findall(hFigSpecr,'Tag','toolbarConvert');...
    findall(hFigSpecr,'Tag','toolbarFootprint');...
    findall(hFigSpecr,'Tag','specr_PushbuttonSelectScan');...
    ];
initEnable = get(hEnableGroup,'Enable');

if isempty(file) | isempty(scan.fidPos)
    set(hToolbarMonitor,'state','off');
    return;
end

icons = load('icons.mat');
if strcmp(get(hToolbarMonitor,'state'),'on');
    set(hToolbarMonitor,'CData',icons.iconmonitor.pause);
else
    set(hToolbarMonitor,'CData',icons.iconmonitor.play);
    if ~isempty(hTimerMonitor)
        if isvalid(hTimerMonitor)
            stop(hTimerMonitor);
        end
        delete(hTimerMonitor);
    end
    return;
end
hTimerMonitor=timer(...
    'Tag','timerMonitor',...
    'StartDelay', 1,...
    'TimerFcn', @timerMonitorFcn,...
    'StartFcn', {@timerMonitorStartFcn,hEnableGroup,initEnable},...
    'StopFcn',  {@timerMonitorStopFcn,hEnableGroup,initEnable},...
    'ExecutionMode','fixedRate',...
    'BusyMode','drop',...
    'TasksToExecute',inf,...
    'Period',settings.monitorPeriod,...
    'userdata',{hEnableGroup,initEnable});
start(hTimerMonitor);

%================================================================
% --- callback function of monitor timer
%================================================================
function timerMonitorFcn(hObject,eventdata)
tStart=tic;
hFigSpecr = findall(0,'Tag','specr_Fig');
hTimerMonitor   = timerfind('Tag','timerMonitor');
hPopupmenuX = findall(hFigSpecr,'Tag','specr_PopupmenuX');
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
settings = getappdata(hFigSpecr,'settings');

file = settings.file;
scan = settings.scan;
% --- check dynamically if there are more scans added to spec file
[fid,message] = fopen(file);
fseek(fid,scan.fidPos(end),'bof');  % move fid to end of the last stored scan
scanline = fgetl(fid);      % skip the currentline (head of the last stored scan);
flagNewScan = 0;
while feof(fid) == 0
    tempfid = ftell(fid);
    scanline = fgetl(fid);
    if length(scanline) >= 3 && strcmp(scanline(1:3),'#S ')
        scan.fidPos = [scan.fidPos,tempfid];
        scan.head{length(scan.fidPos)} = scanline;
        flagNewScan = 1;
    end
end
fclose(fid);
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);
settings = getappdata(hFigSpecr,'settings');
selection = length(scan.fidPos);

if settings.monitorErasemode == 1 | ~isfield(scan,'selection')
    scan.selection = cell(length(selection),1);
    scan.selection{1,1}.time = '';               % initialize scan.selection
    scan.selectionIndex = selection;        % store the selected scan index
    scan.selectionNumber = [];
else
    if scan.selectionIndex(end) ~= selection
        scan.selection(end+1,1) = cell(1,1);
        scan.selectionIndex = [scan.selectionIndex,selection];
    end
end
tmp_selectionNumber = [];

[fid,message] = fopen(file);
try
    fseek(fid,scan.fidPos(selection),'bof');
    % --- read the scan number      PS: scan number sometimes is not the
    % same of scan index.
    scanline = fgetl(fid);
    while ~strcmp(scanline(1:3),'#S ')
        scanline = fgetl(fid);
    end
    space_pos = findstr(scanline,' ');
    tmp_selectionNumber = str2num(scanline(space_pos(1)+1:space_pos(2)-1));
    if isempty(scan.selectionNumber) | scan.selectionNumber(end) ~= tmp_selectionNumber
        scan.selectionNumber = [scan.selectionNumber, tmp_selectionNumber];
    end
    % --- read time of scan #D
    while ~strcmp(scanline(1:2),'#D')
        scanline = fgetl(fid);
    end
    scan.selection{end,1}.time = scanline(4:end);
    % --- read current motor positions
    scan.selection{end,1}.motorPos = [];
    while ~strcmp(scanline(1:2),'#P')
        scanline = fgetl(fid);
    end
    while length(scan.selection{end,1}.motorPos) < length(scan.motorName)
        tmpspace = findstr(scanline,' ');
        scanline = scanline(tmpspace(1):end);
        scan.selection{end,1}.motorPos = [scan.selection{end,1}.motorPos, str2num(scanline)];
        scanline = fgetl(fid);
    end
    % --- read number of scan columns #N ...
    while ~strcmp(scanline(1:2),'#N')
        scanline = fgetl(fid);
    end
    scan.selection{end,1}.length = str2num(scanline(3:end));
    % --- read head of scan columns #L ...
    while ~strcmp(scanline(1:2),'#L')
        scanline = fgetl(fid);
    end
    scanline = scanline(4:end);
    space = findstr(scanline,'  ');
    lengthSpace = length(space);
    for iSpace = lengthSpace:-1:2
        if space(iSpace) == space(iSpace-1)+1
            space(iSpace) = [];
        end
    end
    space = [-1 space length(scanline)+1];
    scan.selection{end,1}.colHead = ...
        cell(1,scan.selection{end,1}.length);
    for iCol = 1:scan.selection{end,1}.length
        colHead = scanline(space(iCol)+2:space(iCol+1)-1);
        while colHead(1) == ' '
            colHead(1) = '';
        end
        scan.selection{end,1}.colHead{iCol} = colHead;
    end

    % --- read scan data
    colData = [];
    str_scanline = num2str(scanline);
    while ~strcmp(str_scanline,'') & ~strcmp(str_scanline,'-1')
        fidPos = ftell(fid);
        scanline = fgetl(fid);
        str_scanline = num2str(scanline);
        while ~strcmp(str_scanline,'') & ~strcmp(str_scanline,'-1')...
                & ~strcmp(str_scanline(1:1),'#')
            fseek(fid,fidPos,'bof');
            colData = [colData fscanf(fid,'%g',[scan.selection{end,1}.length,1])];
            scanline = fgetl(fid);
            fidPos = ftell(fid);
            scanline = fgetl(fid);
            str_scanline = num2str(scanline);
        end
    end
    scan.selection{end,1}.colData = colData';
catch
    return;
end
fclose(fid);
% --- if not the same type scan selected, switch back to erase mode
if length(scan.selection) >= 2
    for iSelection = 2:length(scan.selection)
        if scan.selection{iSelection}.length == scan.selection{iSelection-1}.length
            colHeadCmp = strcmp(scan.selection{iSelection}.colHead,...
                scan.selection{iSelection-1}.colHead);
        else
            colHeadCmp = 0;
        end
        if ~isempty(find(colHeadCmp == 0))
            scan.selection(1:end-1) = '';
            scan.selectionIndex(1:end-1) = '';
            scan.selectionNumber(1:end-1) = '';
            break;
        end
    end
end
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);

if get(hTimerMonitor,'TasksExecuted') == 1 | flagNewScan == 1;
    set(hPopupmenuX,...
        'String',scan.selection{end}.colHead,...
        'Value',1);
    ihklscan = findstr(lower(scan.head{end}),' hklscan');    % check if it is hklscan
    if ihklscan
        a = scan.head{end};
        b = str2num(a(ihklscan+8:end));
        c = find(abs(b([1,3,5])-b([2,4,6]))>eps);
        if c
            set(hPopupmenuX,...
                'String',scan.selection{end}.colHead,...
                'Value',c(1));
        end
    end
    set(hPopupmenuY,...
        'String',scan.selection{end}.colHead,...
        'Value',scan.selection{end}.length);
end
try
    scanplot;
catch
end
tElapsed = ceil(toc(tStart)*1000)/1000;
% --- set checking period
if settings.monitorAutoPeriodMode ~= 1 & tElapsed > settings.monitorPeriod
    settings.monitorAutoPeriodMode = 1;
    settings.monitorAutoPeriod = settings.monitorPeriod;
end
if settings.monitorAutoPeriodMode == 1;
    if tElapsed > settings.monitorAutoPeriod
        settings.monitorAutoPeriod = tElapsed+0.1;
        setappdata(hFigSpecr,'settings',settings);
        userData = get(hTimerMonitor,'userdata');
        hEnableGroup = userData{1};
        initEnable = userData{2};
        stop(hTimerMonitor);
        delete(hTimerMonitor);
        hTimerMonitor=timer(...
            'Tag','timerMonitor',...
            'StartDelay', 1,...
            'TimerFcn', @timerMonitorFcn,...
            'StartFcn', {@timerMonitorStartFcn,hEnableGroup,initEnable},...
            'StopFcn',  {@timerMonitorStopFcn,hEnableGroup,initEnable},...
            'ExecutionMode','fixedRate',...
            'BusyMode','drop',...
            'TasksToExecute',inf,...
            'Period',max(settings.monitorAutoPeriod,settings.monitorPeriod),...
            'userdata',{hEnableGroup,initEnable});
         start(hTimerMonitor);
    end
end


%================================================================
% --- start callback function of monitor timer
%================================================================
function timerMonitorStartFcn(hObject,eventdata,hEnableGroup,initEnable)
hFigSpecr = findall(0,'Tag','specr_Fig');
resettoolbar(hFigSpecr);
set(hEnableGroup,'Enable','off');


%================================================================
% --- stop callback function of monitor timer
%================================================================
function timerMonitorStopFcn(hObject,eventdata,hEnableGroup,initEnable)
for i=1:length(hEnableGroup)
    set(hEnableGroup(i),'Enable',initEnable{i});
end

