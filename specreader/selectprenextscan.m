function selectprenextscan
hTool = gcbo;

hFigSpecr = findall(0,'Tag','specr_Fig');
updatespecscan(hFigSpecr);

settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
if isempty(file) | isempty(scan.fidPos)
    return;
end

% get selection
if ~isfield(scan,'selectionIndex')
    selection = 1;
else
    if strcmpi(get(hTool,'tag'),'toolbarScanPre')
        selection = max(1,scan.selectionIndex(1)-1);
    else
        selection = min(length(scan.fidPos),scan.selectionIndex(end)+1);
    end
end

% read data
scan = readspecscan(file,scan,selection);
% --- post selection scan
postselectscan(hFigSpecr,scan,selection);