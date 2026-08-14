function postselectscan(hFigSpecr,scan,selection)
% --- if not the same type scan selected, error and return; otherwise save
% data
if length(selection) >= 2
    for iSelection = 2:length(selection)
        if scan.selection{iSelection}.length == scan.selection{iSelection-1}.length
            colHeadCmp = strcmp(scan.selection{iSelection}.colHead,...
                scan.selection{iSelection-1}.colHead);
        else
            colHeadCmp = 0;
        end
        if ~isempty(find(colHeadCmp == 0))
            uiwait(msgbox('Multi-selections have to be scans of the same type.',...
                'Select Scan Error','error','modal'));
            return;
        end
    end
end
scan.selectionIndex = selection;        % store the selected scan index (can be different
                                        % from the real scan number is spec file, because of
                                        % the reset of scan number in spec)
settings = getappdata(hFigSpecr,'settings');
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);

% prepare for plotting
hPopupmenuX = findall(hFigSpecr,'Tag','specr_PopupmenuX');
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
hPopupmenuPlotStyle = findall(hFigSpecr,'Tag','specr_PopupmenuPlotStyle');
set(hPopupmenuX,...
    'String',scan.selection{end}.colHead,...
    'Value',1);

ihklscan = findstr(lower(scan.head{scan.selectionIndex(end)}),' hklscan');    % check if it is hklscan
if ihklscan
    a = scan.head{scan.selectionIndex(end)};
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

try
    scanplot;
catch
end
zoom out;