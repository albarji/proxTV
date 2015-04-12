### Example script showing how to perform a 3D Total-Variation filtering with proxTV
import prox_tv as ptv
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import time
import skimage as ski
from skimage import data, io, filters, color, util

# Load color image (3 dimensions: length, width and color)
X = io.imread('colors.png')
X = ski.img_as_float(X)

# Introduce noise
noiseLevel = 0.2
N = util.random_noise(X, mode='gaussian', var=noiseLevel)

# Filter using 3D TV-L1: we only want to penalize X and Y dims (not color)
lam=100./255.;
print('Filtering image...')
start = time.time()
F = ptv.tvgen(N,      [lam, lam],                  [1, 2],                   [1, 1])
#             Image | Penalty in each dimension |  Dimensions to penalize  | Norms to use
end = time.time()
print('Elapsed time ' + str(end-start))

# Any dimension can be penalized under any norm. By also penalizing the color dimension under TV-L2 we get a "decolored" image
lam2=50./255.;
print('Color filtering...')
start = time.time()
F2 = ptv.tvgen(N,      [lam, lam, lam2],            [1, 2, 3],                [1, 1, 2]);
#              Image | Penalty in each dimension |  Dimensions to penalize  | Norms to use 
end = time.time()
print('Elapsed time ' + str(end-start))

# Plot results
plt.subplot(2, 2, 1)
io.imshow(X)
plt.title('Original')

plt.subplot(2, 2, 2)
io.imshow(N)
plt.title('Noisy')

plt.subplot(2, 2, 3)
io.imshow(F)
plt.title('Filtered')

plt.subplot(2, 2, 4)
io.imshow(F2)
plt.title('Color filtered')

plt.show()

