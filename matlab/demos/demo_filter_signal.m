%%% Example script showing how to perform a Total-Variation filtering with proxTV

clear all
close all

%%% TV-L1 filtering

% Generate impulse (blocky) signal
N = 1000;
s = zeros(N,1);
s(N/4:N/2) = 1;
s(N/2:3*N/4) = -1;
s(3*N/4:end-N/8) = 2;

% Introduce noise
n = s + 0.5*randn(size(s));

% Filter using TV-L1
lambda=20;
disp('Filtering signal...');
tic;
f = TV(n,lambda);
toc;

% Plot results
figure();
subplot(3,1,1);
plot(s);
title('Original');
grid();
subplot(3,1,2);
plot(n);
title('Noisy');
grid();
subplot(3,1,3);
plot(f);
title('TVL1-filtered');
grid();

%%% TV-L2 filtering

% Generate sinusoidal signal
N = 1000;
s = sin((1:N)./10) + sin((1:N)./100);

% Introduce noise
n = s + 0.5*randn(size(s));

% Filter using TV-L2
lambda=100;
disp('Filtering signal...');
tic;
f = TV(n,lambda,2);
toc;

% Plot results
figure();
subplot(3,1,1);
plot(s);
title('Original');
grid();
subplot(3,1,2);
plot(n);
title('Noisy');
grid();
subplot(3,1,3);
plot(f);
title('TVL2-filtered');
grid();

%%% Weighted TV-L1 filtering

% Generate impulse (blocky) signal
N = 1000;
s = zeros(N,1);
s(N/4:N/2) = 1;
s(N/2:3*N/4) = -1;
s(3*N/4:end-N/8) = 2;

% Introduce noise
n = s + 0.5*randn(size(s));

% Generate weights
lambda = linspace(0,2,N-1);

% Filter using TV-L1
disp('Filtering signal...');
tic;
f = TV(n,lambda);
toc;

% Plot results
figure();
subplot(4,1,1);
plot(s);
title('Original');
grid();
subplot(4,1,2);
plot(n);
title('Noisy');
grid();
subplot(4,1,3);
plot(f);
title('Weighted TVL1-filtered');
grid();
subplot(4,1,4);
area(lambda);
title('Weights');
grid();
