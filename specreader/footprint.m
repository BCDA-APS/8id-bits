function footprint(varargin)
% FOOTPRINTGUI Called by specr and trueref.
%
% Copyright 2004, Zhang Jiang

hFigSpecr = findall(0,'Tag','specr_Fig');
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

% --- get wavelength and footprint angle from main figure
settings = getappdata(hFigSpecr,'settings');
footprintAngle = settings.footprintAngle;
wavelength = settings.wavelength;

% --- footprint correction
for iLine = 1:length(hLine)
    xdata = get(hLine(iLine),'XDATA');
    ydata = get(hLine(iLine),'YDATA');
    ydataError = getappdata(hLine(iLine),'ydataError');
    rawData = [xdata' ydata' ydataError];
    footprintData = footprintFcn(rawData,footprintAngle,wavelength);
    xdata = footprintData(:,1)';
    ydata = footprintData(:,2)';
    ydataError = footprintData(:,3);
    set(hLine(iLine),'XDATA',xdata);
    set(hLine(iLine),'YDATA',ydata);
    setappdata(hLine(iLine),'ydataError',ydataError);
end

% ======================================================================
% --- calculated footprint correted data
% ======================================================================
function result = footprintFcn(varargin)
% FOOTPRINTFCN Footprint correction to reflectivity
%   RESULT = FOOTPRINTFCN(RAWDATA,FOOTANGLE,WAVELENGTH)
%
%   Format of input:
%       RAWDATA: Mx2 or Mx3 matrix. 1st column qz (Unit: A^-1), 2nd 
%           intensity, 3rd absolute error of intensity if exists.
%       FOOTANGLE: footprint angle (Unit: degree)
%       WAVELENGTH: beam wavelength (Unit: angstrom)
%
%   Format of output:
%       RESULT: Mx2 or Mx3
%
% Copyright 2004 Zhang Jiang
% $Revision: 1.0 $  $Date: 2004/11/19 11:35:28 $

rawData = varargin{1};      % M x 2(3)with 1st col qz, 2nd int, 3rd abs error
footAngle = varargin{2};    % in degree
wavelength = varargin{3}; % Angstrom

footQz = abs(4*pi/wavelength*sin(footAngle*pi/180));
correctIndex = find(abs(rawData(:,1))<footQz);

% --- correct intensity
result = rawData;
warning off MATLAB:divideByZero;
if ~isempty(correctIndex)
    result(correctIndex,2) = ...
        rawData(correctIndex,2)*footQz./abs(rawData(correctIndex,1));

end

% --- correction to the absolute error
if size(rawData,2) == 3     
    result(:,3) = result(:,2).*(rawData(:,3)./rawData(:,2));
end
warning on MATLAB:divideByZero;