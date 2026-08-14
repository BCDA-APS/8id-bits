function meshscanplot
hFigSpecr = findall(0,'Tag','specr_Fig');
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
if isempty(file) | isempty(scan) | ~isfield(scan,'selectionIndex') | isempty(scan.selectionIndex)
    return;
end

% --- get z data coloum
hPopupmenuY = findall(hFigSpecr,'Tag','specr_PopupmenuY');
popupmenuYString = get(hPopupmenuY,'String');
popupmenuYValue = get(hPopupmenuY,'Value');

% --- load data and plot
[~,filename,~] = fileparts(file);
for iscan = 1:length(scan.selectionIndex)
    scanParts = regexp(scan.head{scan.selectionIndex(iscan)}, '\ *', 'split');
    if ~strcmpi(scanParts{3},'mesh')
        continue;
    end
    xdata = scan.selection{iscan}.colData(:,1);
    ydata = scan.selection{iscan}.colData(:,2);
    zdata = scan.selection{iscan}.colData(:,popupmenuYValue);
    ptsX = str2double(scanParts{7})+1;
    ptsY = floor(length(xdata)/ptsX);
    nOfPts = ptsX*ptsY;
    if nOfPts == 0
        continue;
    end
    meshdata.X=reshape(xdata(1:nOfPts), ptsX, ptsY)';
    meshdata.Y=reshape(ydata(1:nOfPts), ptsX, ptsY)';
    meshdata.Z=reshape(zdata(1:nOfPts), ptsX, ptsY)';
    title_str = [filename,' - Scan #',num2str(scan.selectionNumber(iscan))];
    figure('Name',title_str);
    switch settings.meshScanPlotType
        case 1
            contour(meshdata.X,meshdata.Y,meshdata.Z);
            colormap;
        case 2
            contourf(meshdata.X,meshdata.Y,meshdata.Z);
            colormap(flipud(colormap('gray')));
        case 3
            mesh(meshdata.X,meshdata.Y,meshdata.Z);
            colormap;
        case 4
            surf(meshdata.X,meshdata.Y,meshdata.Z);
            colormap;           
    end
    xlabel(titlestr(scan.selection{iscan}.colHead{1}));
    ylabel(titlestr(scan.selection{iscan}.colHead{2}));
    %set(gca,'String',popupmenuYString(popupmenuYValue));    
    title(titlestr(title_str));
    hC = colorbar;
    ylabel(hC,titlestr(popupmenuYString{popupmenuYValue}));
    if settings.meshScanExport2WS == 1
        varname = ['specrMeshScan',num2str(scan.selectionNumber(iscan))];
        assignin('base',varname,meshdata);
    end
end 
return;