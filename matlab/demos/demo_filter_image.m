%%% Example script showing how to perform a 2D Total-Variation filtering with proxTV

clear all
close all

% Load image
X = rgb2gray(imread('colors.png'));
X=double(X)/255;

% Introduce noise
noiseLevel = 0.01;
N = double(imnoise(X,'speckle'));

% Filter using 2D TV-L1
lambda=0.15;
disp('Filtering image...');
tic;
F = TV(N,lambda);
toc;

% Plot results
figure();
colormap gray
subplot(1,3,1);
imagesc(X);
title('Original');
subplot(1,3,2);
imagesc(N);
title('Noisy');
subplot(1,3,3);
imagesc(F);
title('Filtered');


