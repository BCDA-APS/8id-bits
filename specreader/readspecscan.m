function scan = readspecscan(file,scan,selection)
scan.selection = cell(length(selection),1);
[fid,message] = fopen(file);
scan.selectionNumber = [];
for iSelection = 1:length(selection)
    fseek(fid,scan.fidPos(selection(iSelection)),'bof');
    % --- read the scan number      PS: scan number sometimes is not the
    % same of scan index.
    scanline = fgetl(fid);
    while ~strcmp(scanline(1:3),'#S ')
        scanline = fgetl(fid);
    end
    space_pos = findstr(scanline,' ');
    scan.selectionNumber(iSelection) = str2num(scanline(space_pos(1)+1:space_pos(2)-1));
    % --- read time of scan #D
    while ~strcmp(scanline(1:2),'#D')
        scanline = fgetl(fid);
    end
    scan.selection{iSelection}.time = scanline(4:end);        
    % --- read current motor positions
    scan.selection{iSelection}.motorPos = [];
    while ~strcmp(scanline(1:2),'#P')
        scanline = fgetl(fid);
    end
    while length(scan.selection{iSelection}.motorPos) < length(scan.motorName)
        tmpspace = findstr(scanline,' ');
        scanline = scanline(tmpspace(1):end);
        scan.selection{iSelection}.motorPos = [scan.selection{iSelection}.motorPos, str2num(scanline)];
        scanline = fgetl(fid);
    end
    % --- read number of scan columns #N ...
    while ~strcmp(scanline(1:2),'#N')   
        scanline = fgetl(fid);
    end
    scan.selection{iSelection}.length = str2num(scanline(3:end));
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
    scan.selection{iSelection}.colHead = ...
        cell(1,scan.selection{iSelection}.length);
    for iCol = 1:scan.selection{iSelection}.length
        colHead = scanline(space(iCol)+2:space(iCol+1)-1);
        while colHead(1) == ' '
            colHead(1) = '';
        end
        scan.selection{iSelection}.colHead{iCol} = colHead;
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
            colData = [colData fscanf(fid,'%g',[scan.selection{iSelection}.length,1])];
            scanline = fgetl(fid);
            fidPos = ftell(fid);
            scanline = fgetl(fid);
            str_scanline = num2str(scanline);
        end
    end
    scan.selection{iSelection}.colData = colData';
end
fclose(fid);