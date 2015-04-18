### Example script showing how to perform a 2D Total-Variation filtering with proxTV
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

# Filter using 2D TV-L1
lam=0.15;
print('Filtering image with 2D TV-L1...')
start = time.time()
F = ptv.tv1_2d(N, lam)
end = time.time()
print('Elapsed time ' + str(end-start))

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

plt.show()

