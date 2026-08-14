function updatespecscan(hFigSpecr)
settings = getappdata(hFigSpecr,'settings');
file = settings.file;
scan = settings.scan;
if isempty(file) | isempty(scan.fidPos)
    return;
end

% check dynamically if there are more scans added to spec file
[fid,message] = fopen(file);
fseek(fid,scan.fidPos(end),'bof');  % move fid to end of the last stored scan
scanline = fgetl(fid);      % skip the currentline (head of the last stored scan);
while feof(fid) == 0
    tempfid = ftell(fid);
    scanline = fgetl(fid);
    if length(scanline) >= 3 && strcmp(scanline(1:3),'#S ')
        scan.fidPos = [scan.fidPos,tempfid];
        scan.head{length(scan.fidPos)} = scanline;
    end
end
fclose(fid);
settings.scan = scan;
setappdata(hFigSpecr,'settings',settings);