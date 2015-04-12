### Example script showing how to perform a weighted 2D Total-Variation filtering with proxTV
import prox_tv as ptv
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import time
import skimage as ski
from skimage import data, io, filters, color, util

# Load image
X = io.imread('colors.png')
X = ski.img_as_float(X)
X = color.rgb2gray(X)

# Introduce noise
noiseLevel = 0.01
N = util.random_noise(X, mode='speckle', var=noiseLevel)

# Gradient in columns
W1 = 0.01 * np.cumsum(np.ones((X.shape[0]-1, X.shape[1])), 1)
W2 = 0.01 * np.ones((X.shape[0], X.shape[1]-1))
print('Solving 2D weighted TV...' )
start = time.time()
FW = ptv.tv1w_2d(N, W1, W2)
end = time.time()
print('Elapsed time ' + str(end-start))

plt.subplot(3, 4, 1); io.imshow(W1); plt.title('Weights along columns');
plt.subplot(3, 4, 5); io.imshow(W2); plt.title('Weights along rows');
plt.subplot(3, 4, 9); io.imshow(FW); plt.title('Filter result');

# Gradient in rows
W1 = 0.01 * np.ones((X.shape[0]-1, X.shape[1]))
W2 = 0.01 * np.cumsum(np.ones((X.shape[0], X.shape[1]-1)), 0)
print('Solving 2D weighted TV...')
start = time.time()
FW = ptv.tv1w_2d(N, W1, W2)
end = time.time()
print('Elapsed time ' + str(end-start))

plt.subplot(3, 4, 2); io.imshow(W1); plt.title('Weights along columns');
plt.subplot(3, 4, 6); io.imshow(W2); plt.title('Weights along rows');
plt.subplot(3, 4, 10); io.imshow(FW); plt.title('Filter result');

# Gradient in columns and rows
W1 = 0.004 * np.cumsum(np.ones((X.shape[0]-1, X.shape[1])), 1)
W2 = 0.004 * np.cumsum(np.ones((X.shape[0], X.shape[1]-1)), 0)
print('Solving 2D weighted TV...' )
start = time.time()
FW = ptv.tv1w_2d(N, W1, W2)
end = time.time()
print('Elapsed time ' + str(end-start))

plt.subplot(3, 4, 3); io.imshow(W1); plt.title('Weights along columns');
plt.subplot(3, 4, 7); io.imshow(W2); plt.title('Weights along rows');
plt.subplot(3, 4, 11); io.imshow(FW); plt.title('Filter result');

# Grid regions
W1=cos(0.05*np.cumsum(np.ones((X.shape[0]-1, X.shape[1])), 1))
W1[W1 < 0] = 0
W2=sin(0.05*np.cumsum(np.ones((X.shape[0], X.shape[1]-1)), 0))
W2[W2 < 0] = 0;
print('Solving 2D weighted TV...')
start = time.time()
FW = ptv.tv1w_2d(N, W1, W2)
end = time.time()
print('Elapsed time ' + str(end-start))

plt.subplot(3, 4, 4); io.imshow(W1); plt.title('Weights along columns');
plt.subplot(3, 4, 8); io.imshow(W2); plt.title('Weights along rows');
plt.subplot(3, 4, 12); io.imshow(FW); plt.title('Filter result');

plt.show()

