%% Example of ECG denoising using CNC FLSA
%
% Ankit Parekh, NYU School of Engineering. 
% Ref.: Convex fused lasso denoising with non-convex regularization and its
%       use for pulse detection.
%       Ankit Parekh and Ivan W. Selesnick. 2015
% 

%% Initialize
clear, close all; clc
printme = @(x) print('-dpdf',x);
rmse = @(y,x) sqrt( sum( (y(:)-x(:)).^2) / numel(y) );
%% Generate synthetic ECG signal

ecg = load('ecgSignal.txt');
fs = 256;                                                                   % Sampling frequency
N = length(ecg);
n = 0:N-1;

rng('default')
sigma = 0.4;                                                                % sigma : noise standard deviation
noise = sigma * randn(N, 1);                                                % noise : white Gaussian noise
data = ecg + noise;                                                         % data : noisy ECG
%% Run fused lasso for pulse detection

lam0 = 0.6;
lam1 = 0.9;
a0 = 0.9/lam0;
a1 = (1-a0*lam0) / (4*lam1);
Nit = 50;
pen = 'atan';

xL1 = soft(tvd(data,N,lam0),lam1);                                          % L1 FLSA
[x, cost] = CNC_FLSA(data, lam0, lam1, a0, a1, Nit, pen);                       % CNC FLSA
%% Plot the denoised ECG signals

start = 0;
finish = 30;
mid = (start+finish)/2;
figure(1), clf
gap = 4;
plot(n/fs,ecg,'k', n/fs, data-gap, 'k', n/fs, x-2*gap, 'k',...
                                                    n/fs, xL1-3*gap, 'k')
box off
ylim([-13 2]); xlim([start finish])
xlabel('Time (s)')
text(mid, 2, 'ECG signal','horizontalalignment','center')
text(mid, -1.5, 'Noisy ECG signal','horizontalalignment','center')
text(mid, -6, 'Denoised ECG using CNC fused lasso',...
    'horizontalalignment','center')
text(mid, -10, 'Denoised ECG using L1 fused lasso',...
    'horizontalalignment','center')

set(gca,'YTick', [-13 -12 -11 -9 -8 -7 -5 -4 -3 -1 0 1],...
    'YTickLabel',[-1 0 1]);
set(gca,'XTick', start:3:finish, 'XTickLabel', 0:3:finish-start)
printme('ECG_Signal_Denoising')
