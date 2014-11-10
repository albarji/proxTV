%%% Example script showing how to perform a 3D Total-Variation filtering with proxTV

clear all
close all

% Load color image (3 dimensions: length, width and color)
X = imread('colors.png');

% Introduce noise
noiseLevel = 0.2;
N = double(imnoise(X,'gaussian',0,noiseLevel));

% Filter using 3D TV-L1: for this one needs to invoke prox_TVgen, since we only want to penalize X and Y dims (not color)
lambda=100;
disp('Filtering image...');
tic;
F = TVgen(N,      [lambda lambda],             [1 2],                    [1 1]);
%              Image | Penalty in each dimension |  Dimensions to penalize  | Norms to use
toc;

% Any dimension can be penalized under any norm. By also penalizing the color dimension under TV-L2 we get a "decolored" image
lambda2=50;
disp('Color filtering...');
tic;
F2 = TVgen(N,      [lambda lambda lambda2],     [1 2 3],                  [1 1 2]);
%               Image | Penalty in each dimension |  Dimensions to penalize  | Norms to use 
toc;

% Plot results
figure();
subplot(2,2,1);
imshow(X);
title('Original');
subplot(2,2,2);
imshow(uint8(N));
title('Noisy');
subplot(2,2,3);
imshow(uint8(F));
title('Filtered');
subplot(2,2,4);
imshow(uint8(F2));
title('Color filtered');


