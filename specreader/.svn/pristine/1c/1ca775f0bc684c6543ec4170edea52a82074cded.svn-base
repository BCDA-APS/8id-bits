function normData = norm2maxi(varargin)
% NORM2MAXI Normalize reflectivity to maximum intensity
%   NORMDATA = NORM2MAXI(RAWDATA,ARGUMENT)
%
%   Format of input:
%       RAWDATA: Mx2 or Mx3 matrix. 1st column qz (Unit: A^-1), 2nd 
%           intensity, 3rd absolute error of intensity if exists.
%       ARGUMENT: can be the intensity to be normalized to, or blank. If
%           there is no argument, then normalize to maximum intensity.
%
%   Format of output:
%       NORMDATA: Mx2 or Mx3
%
% Copyright 2004 Zhang Jiang
% $Revision: 1.0 $  $Date: 2004/11/19 17:20:54 $

rawData = varargin{1};      % M x 2(3)with 1st col qz, 2nd int, 3rd abs error

% --- determine normalizing to what intensity
if nargin == 1
    maxI = max(rawData(:,2));
else
    maxI = varargin{2};
end

normData = rawData;
normData(:,2:3) = rawData(:,2:3)/maxI;
    
    
