### Example script showing how to run proxTV solvers in multiple parallel threads
import prox_tv as ptv
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import time
import skimage as ski
from skimage import data, io, filters, color, util

# WRITE HERE YOUR NUMBER OF THREADS
threads = 2

# Load image
X = io.imread('QRbig.png')
X = ski.img_as_float(X)
X = color.rgb2gray(X)

# Introduce noise
noiseLevel = 0.2
N = util.random_noise(X, mode='gaussian', var=noiseLevel)

# Filter using 1 thread
lam=50./255.;
print('Filtering image with 1 thread...');
start = time.time()
F = ptv.tv1_2d(N, lam)
end = time.time()
print 'Elapsed time ' + str(end-start)

# Now filter using several threads
print('Filtering image with ' + str(threads) + ' threads...');
start = time.time()
F = ptv.tv1_2d(N, lam, n_threads=threads)
end = time.time()
print 'Elapsed time ' + str(end-start)

# Plot results
plt.subplot(1, 3, 1)
io.imshow(X)
plt.xlabel('Original')

plt.subplot(1, 3, 2)
io.imshow(N)
plt.title('2D TVL1 filtering')
plt.xlabel('Noisy')

plt.subplot(1, 3, 3)
io.imshow(F)
plt.xlabel('Filtered')

show()

