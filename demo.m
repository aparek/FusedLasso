%% Example of 1D signal denoising using CNC FLSA
%
% Ankit Parekh, NYU School of Engineering. 
% Ref.: Convex fused lasso denoising with non-convex regularization and its
%       use for pulse detection.
%       Ankit Parekh and Ivan W. Selesnick. 2015
% 

%% Initialize
clear, clc, close all;
printme = @(x) print('-dpdf',x);
rmse = @(y,x) sqrt( sum( (y(:)-x(:)).^2) / numel(y) );
%% Generate synthetic test signals

N = 300;
n = 0:N-1;
sigma = 1.5;
rng('default')

s = zeros(N,1);
s(50:55) = -2.2; s(100:105) = -2; s(180:210) = 3; 
s(120:130) = 2.25;
y = s + sigma*randn(size(s));

figure(1)
clf
subplot(2,1,1)
plot(s,'k')
title('Clean Test Signal')

subplot(2,1,2)
plot(y,'k')
title(sprintf('Noisy Test Signal(\\sigma = %1.1f). RMSE = %1.4f',...
    sigma, rmse(s,y)))
xlabel('Time (n)')
printme('Test_Signals')

%% Denoise using CNC and L1 FLSA and plot the cost function history
lam0 = 0.8*sigma;
lam1 = 0.15 * sqrt(N) *sigma;
a0 = 0.9 / lam0;
a1 = (1-a0 * lam0) / (4 * lam1);
Nit = 20;

[xE,cost] = CNC_FLSA(y,lam0,lam1,a0,a1,Nit,'atan');                             %CNC FLSA
xL1 = soft(tvd(y,N,lam1),lam0);                                             %L1 FLSA
xE = xE(:);

figure(2), clf
plot(cost,'.-k')
title('Cost function history for CNC FLSA algorithm using Majorization-Minimization')
xlabel('Iteration number')
ylabel('Value of objective function')
printme('Cost_function_history')

%% Plot figures
figure(3)
clf

subplot(4,1,1)
plot(s,'k')
title('Synthetic Pulse Signal')


subplot(4,1,2)
plot(y,'k')
title('Noisy Pulse Signal (\sigma = 1.5)')

subplot(4,1,3)
plot(n, xE,'k'); hold on; plot(n,xL1,'b');
title('Estimated Pulse Signal')
legend('CNC','L1')

subplot(4,1,4)
plot(n, s-xE,'k'); hold on; plot(n,s-xL1','b');
title('Denoising Error')
legend('CNC','L1')
xlabel('Time (n)')

printme('Denoising_Results')

