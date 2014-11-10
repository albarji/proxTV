%%% Example script showing how to run proxTV solvers in multiple parallel threads

clear all
close all

% WRITE HERE YOUR NUMBER OF THREADS
threads = 2;

% Load image
X = rgb2gray(imread('QRbig.png'));

% Introduce noise
noiseLevel = 0.2;
N = double(imnoise(X,'gaussian',0,noiseLevel));

% Filter using 1 thread
lambda=50;
disp('Filtering image...');
tic;
F = TV(N,lambda);
toc;

% Now filter using several threads
lambda=50;
disp('Filtering image...');
tic;
F = TV(N,lambda,1,threads);
toc;

% Plot results
figure();
subplot(1,3,1);
imshow(X);
title('Original');
subplot(1,3,2);
imshow(uint8(N));
title('Noisy');
subplot(1,3,3);
imshow(uint8(F));
title('Filtered');

