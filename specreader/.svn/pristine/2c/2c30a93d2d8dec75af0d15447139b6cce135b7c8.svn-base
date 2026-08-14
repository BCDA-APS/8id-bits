function varargout = xrr_fft(varargin)
% XRR_FFT Calculate and plot the Fast Fourier Transform (FFT) of 
%   x-ray reflectivity to estimate film thickness.
%   XRR_FFT(DATA) or FFTDATA = XRR_FFT(DATA)
%
%   Format of input:
%       DATA: Mx2 or Mx3 matrix. 1st column qz (Unit: A^-1), 2nd 
%           intensity, 3rd absolute error of intensity if exists.
%
%   Format of output:
%       FFTDATA: Mx2. 1st column is the depth (A); 2nd column is the FFT.
%
% Copyright 2011 Zhang Jiang
% $Revision: 1.0 $  $Date: 2011/11/14 $

if nargin ~= 1 
    error('Invalid inuput argument.');
end

data = varargin{1};
data = sortrows(data,1);
qz = data(:,1);
I = data(:,2);
I = I/max(I);

I_norm = log10(I.*qz.^4);
I_norm = I_norm - mean(I_norm);

L=10*length(qz);
NFFT = 2^nextpow2(L);
Y = fft(I_norm,NFFT)/L;
f = 2*pi/mean(diff(qz))/2*linspace(0,1,NFFT/2+1); f = f';
depth = 2*abs(Y(1:NFFT/2+1));

figure
plot(f,depth);
box on;
xlabel('Depth (A)');
ylabel('|FFT|');

if nargout == 1
    varargout{1} = [f,depth];
end



