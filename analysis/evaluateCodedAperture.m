
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
