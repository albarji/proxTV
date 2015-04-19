%%% Example script showing how to perform a weighted 2D Total-Variation filtering with proxTV

clear all
close all

% Prepare plot
figure()
colormap gray

% Load image
X = rgb2gray(imread('colors.png'));
X=double(X)/255;

% Introduce noise
noiseLevel = 0.01;
N = double(imnoise(X,'speckle'));

% Gradient in columns
W1=0.01*cumsum(ones(size(X,1)-1, size(X,2)),2);
W2=0.01*ones(size(X,1), size(X,2)-1);
disp('Solving 2D weighted TV...');
tic;
FW = TV(N,{W1,W2});
toc;
subplot(3,4,1); imagesc(W1); title('Weights along columns');
subplot(3,4,5); imagesc(W2); title('Weights along rows');
subplot(3,4,9); imagesc(FW); title('Filter result');

% Gradient in rows
W1=0.01*ones(size(X,1)-1, size(X,2));
W2=0.01*cumsum(ones(size(X,1), size(X,2)-1),1);
disp('Solving 2D weighted TV...');
tic;
FW = TV(N,{W1,W2});
toc;
subplot(3,4,2); imagesc(W1); title('Weights along columns');
subplot(3,4,6); imagesc(W2); title('Weights along rows');
subplot(3,4,10); imagesc(FW); title('Filter result');

% Gradient in columnas and rows
W1=0.004*cumsum(ones(size(X,1)-1, size(X,2)),2);
W2=0.004*cumsum(ones(size(X,1), size(X,2)-1),1);
disp('Solving 2D weighted TV...');
tic;
FW = TV(N,{W1,W2});
toc;
subplot(3,4,3); imagesc(W1); title('Weights along columns');
subplot(3,4,7); imagesc(W2); title('Weights along rows');
subplot(3,4,11); imagesc(FW); title('Filter result')

% Grid regions
W1=max(0,cos(0.05*cumsum(ones(size(X,1)-1, size(X,2)),2)));
W2=max(0,sin(0.05*cumsum(ones(size(X,1), size(X,2)-1),1)));
disp('Solving 2D weighted TV...');
tic;
FW = TV(N,{W1,W2});
toc;
subplot(3,4,4); imagesc(W1); title('Weights along columns');
subplot(3,4,8); imagesc(W2); title('Weights along rows');
subplot(3,4,12); imagesc(FW); title('Filter result')

