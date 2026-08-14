function sresult = mergescan(varargin)
% MERGESCAN Merge and plot scans from the same spec file.
%
%   Usage:
%       M = MERGESCAN(FILENAME,SCANS,BACKGROUND) 
%       M = MERGESCAN(FILENAME,SCANS,BACKGROUND,DETECTOR)
%       M = MERGESCAN(FILENAME,SCANS,BACKGROUND,PLOTSWITCH)
%       M = MERGESCAN(FILENAME,SCANS,BACKGROUND,DETECTOR,PLOTSWITCH)
%
%   Input Argument:
%       FILENAME    : Spec file name string.
%       SCANS       : List of scans to be merged, e.g., [17,18,19].
%       BACKGROUND  : List of background of corresponding backgournd, e.g., 
%                       [110,110,110]. Enter [0,0,0] for zero background.
%       DETECTOR    : Detecor string; default value is the last column
%       PLOTSWITCH  : {'on'} or {'off'}; default is {'off'}. When it is 
%                       {'on'}, plot style can be specified in form of 
%                       {'on',LineSpec}, where LineSpec is defined in 
%                       matlab buildin function ERRORBAR. For more
%                       information on LineSpec, refer to ERRORBAR. The plot
%                       is set to logy scale automatically.
%
%   Out Argument:
%       M : three column matrix, the 3rd col of which is the absolute
%           statistical error.
%
%   Example:
%       m = mergescan('ba5103ps',[17,18,19],[370,370,370])
%       m = mergescan('ba5103ps',[17,18,19],[370,370,370],'pind3')
%       m = mergescan('ba5103ps',[17,18,19],[370,370,370],{'on'})
%       m = mergescan('ba5103ps',[17,18,19],[370,370,370],{'on','ro-'})
%       m = mergescan('ba5103ps',[17,18,19],[370,370,370],'pind3',{'on','k^:'})
%
% Copyright 2006 Zhang Jiang
% $Revision: 1.0 $  $Date: 2006/01/24 $

if nargin < 3 | nargin > 5
    error('Invalid input argument.');
    return;
end
f           = varargin{1};
scan        = varargin{2};
bkg         = varargin{3};
if ~ischar(f)
    error('Invalid spec filename.');
    return;
end
if ~isnumeric(scan) | min(size(scan)) ~= 1
    error('Invalid scan numbers.');
    return;
end
if ~isnumeric(bkg) | numel(scan)-numel(bkg) ~= 0
    error('Invalid background or background dimension differs from scan dimension.');
    return;
end
switch nargin
    case 3
        detector = '-1';
        plotswitch = {'off'};
    case 4
        if ischar(varargin{4})
            detector = varargin{4};
            plotswitch = {'off'};
        elseif iscell(varargin{4})
            detector = '-1';
            plotswitch = varargin{4};
        else
            error('Invalid input argument.')'
            return;
        end
    case 5
        detector      = varargin{4};
        plotswitch    = varargin{5};
end
if length(plotswitch) == 1
    if ~ischar(plotswitch{1}) | (~strcmpi(plotswitch{1},'off') & ~strcmpi(plotswitch{1},'on'))
        error('Invalid plotswitch.');
        return;
    end
elseif length(plotswitch) == 2
    if ~ischar(plotswitch{1}) | ~strcmpi(plotswitch{1},'on')...
            | ~ischar(plotswitch{2})
        error('Invalid plotswitch.');
        return;
    end
else
    error('Two many arguments in plotswitch.');
    return;
end

% merge settings
merge.mode = 1;             % intersnity based
merge.interpMethod = 1;     % interplation method

% --- load data
nOfScan     = length(scan);
scanData    = cell(1,nOfScan);
for iOfScan = 1:nOfScan
    tempScan = rdspec(f,scan(iOfScan));
    if strcmp('-1',detector);
        headIndex = length(tempScan.head);
    else
        headIndex = 0;
        for iOfHead = 1:length(tempScan.head)
            if strcmpi(detector,tempScan.head{iOfHead})
                headIndex = iOfHead;
            end
        end
        if headIndex == 0
            error('Invalid detector name.');
            return;
        end
    end
    xdata = tempScan.data(:,1);
    ydata = tempScan.data(:,headIndex);
    ydataRelErr = 1./sqrt(ydata);
    ydata = ydata - bkg(iOfScan);
    ydataErr = ydata.*ydataRelErr;
    negIndex = find(ydata<=0);
    xdata(negIndex)         = [];
    ydata(negIndex)         = [];
    ydataErr(negIndex)      = [];    
    scanData{iOfScan} = [xdata,ydata,ydataErr];
    clear tempScan headIndex iOfHead
end

% --- mergescans
sresult = mrgker(scanData,merge);

% --- plot
if strcmpi(plotswitch{1},'off')
    return;
elseif strcmpi(plotswitch{1},'on') & length(plotswitch)==1
    errorbar(sresult(:,1),sresult(:,2),sresult(:,3));
elseif strcmpi(plotswitch{1},'on') & length(plotswitch)==2
    errorbar(sresult(:,1),sresult(:,2),sresult(:,3),plotswitch{2});
else
    return;
end
figure(gcf);
set(gca,'yscale','log');
box on;

% -- EOF