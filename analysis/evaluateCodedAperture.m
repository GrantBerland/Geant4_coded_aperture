
load('./CA_files/mask.mat', 'mask')
load('./CA_files/decoder.mat', 'decoder')

%%
% $SPSF = A * G$
SPSF = conv2(mask, decoder);

%% 
% $SPSF = A * G = \bar{W} W$
% 
% In some cases, $\bar{W} = W^T = W^{-1}$

W = sqrtm(SPSF);

%%
% $MTF = F_d(A * G) = F_d(A)F_d(G)$
MTF = fft2(SPSF);

figure();
imagesc(abs(fftshift(MTF)));
colorbar(); title('Fourier Transform of Point Spread Function');

%%
% $\sigma^2 = (S_{tot}+B_{tot})\tau^2 (Trace(WW^T))^{-1}$
%
tr_inv = 1/trace(W*W');

tau = 145/289; % M/N, ratio of open elements to total elements, N = 2M - 1
S = 100;       % Signal counts on detector
B = 1000;      % Background counts on detector

sigma2 = (S + B) * tau^2 * tr_inv;

SNR = S/sqrt(sigma2);

%%

pixelSize    = 2.2;   % mm
pixelPitch   = 2.46;  % mm
detectorSize = 39.12; % mm
nPixels      = 256;

a = pixelSize * sqrt(nPixels);
b = pixelSize * sqrt(nPixels);

c = detectorSize;
d = detectorSize;

% Fill factor
FF = (a*b)/(c*d);

% Comb function

x = linspace(0, detectorSize, detectorSize/(pixelPitch-1));
diraccomb = zeros(size(x));
diraccomb(1) = 1;
diraccomb = repmat(diraccomb,1,10);

twoDcomb  = diraccomb'*diraccomb;

% Rect function

r1 = rectangularPulse(-pixelSize/2, pixelSize/2, linspace(-pixelPitch/2, pixelPitch/2));
pixelRow = repmat(r1, 1, 16);

twoDrect = pixelRow'*pixelRow;

a = repmat(rectangularPulse(-pixelSize/2, pixelSize/2, linspace(-detectorSize/2, detectorSize/2)), 1, 16);
twoDrect2 = a'*a;

sampling_fnc = conv2(twoDcomb, twoDrect, 'same') * twoDrect2;

sensorMTF = fft2(sampling_fnc);

figure()
contour(real(fftshift(sensorMTF)));
