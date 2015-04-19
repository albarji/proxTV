%%% Example script showing how to run proxTV solvers in multiple parallel threads

clear all
close all

% WRITE HERE YOUR NUMBER OF THREADS TO TEST
THREADS = [1 2 3 4 5 6 7 8];

% Load image
X = rgb2gray(imread('QRbig.png'));

% Introduce noise
noiseLevel = 0.2;
N = double(imnoise(X,'gaussian',0,noiseLevel));

% Iterate over number of threads
times = zeros(size(THREADS));
idx = 1;
for nthreads = THREADS
    % Filter image
    lambda=50;
    disp(['Filtering image with ', num2str(nthreads), ' threads ...']);
    tic;
    F = TV(N, lambda, 1, nthreads);
    times(idx) = toc;
    disp([num2str(toc), ' seconds']);
    idx = idx + 1;
end

% Plot filtering results
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

% Plot measured times
figure();
bar(THREADS, times);
xlabel('Threads','FontSize',20);
ylabel('Time (s)','FontSize',20);
title('Filtering times for increasing threads','FontSize',25);
set(gca,'FontSize',17);
