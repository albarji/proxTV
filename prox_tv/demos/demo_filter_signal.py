### Example script showing how to perform a Total-Variation filtering with proxTV
import prox_tv as ptv
import numpy as np
import matplotlib.pyplot as plt
import time

def _blockysignal():
    """Generates a blocky signal for the demo"""
    N = 1000
    s = np.zeros((N,1))
    s[int(N/4):int(N/2)] = 1
    s[int(N/2):int(3*N/4)] = -1
    s[int(3*N/4):int(-N/8)] = 2
    return s

### TV-L1 filtering

# Generate impulse (blocky) signal
s = _blockysignal()

# Introduce noise
n = s + 0.5*np.random.rand(*np.shape(s))

# Filter using TV-L1
lam=20
print('Filtering signal with TV-L1...')
start = time.time()
f = ptv.tv1_1d(n,lam)
end = time.time()
print('Elapsed time ' + str(end-start))

# Plot results
plt.subplot(3, 1, 1)
plt.title('TVL1 filtering')
plt.plot(s)
plt.ylabel('Original')
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(n)
plt.ylabel('Noisy')
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(f)
plt.ylabel('Filtered')
plt.grid(True)

plt.show()

### TV-L2 filtering

# Generate sinusoidal signal
N = 1000
s = np.sin(np.arange(1,N+1)/10.0) + np.sin(np.arange(1,N+1)/100.0)

# Introduce noise
n = s + 0.5*np.random.randn(*np.shape(s))

# Filter using TV-L2
lam=100;
print('Filtering signal with TV-L2...')
start = time.time()
f = ptv.tv2_1d(n,lam);
end = time.time()
print('Elapsed time ' + str(end-start))

# Plot results
plt.subplot(3, 1, 1)
plt.title('TVL2 filtering')
plt.plot(s)
plt.ylabel('Original')
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(n)
plt.ylabel('Noisy')
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(f)
plt.ylabel('Filtered')
plt.grid(True)

plt.show()

### Weighted TV-L1 filtering

# Generate impulse (blocky) signal
s = _blockysignal()

# Introduce noise
n = s + 0.5*np.random.randn(*np.shape(s))

# Generate weights
lam = np.linspace(0,2,N-1)

# Filter using weighted TV-L1
print('Filtering signal with weighted TV-L1...')
start = time.time()
f = ptv.tv1w_1d(n, lam)
end = time.time()
print('Elapsed time ' + str(end-start))

# Plot results
plt.subplot(4, 1, 1)
plt.title('Weighted TVL1 filtering')
plt.plot(s)
plt.ylabel('Original')
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(n)
plt.ylabel('Noisy')
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(f)
plt.ylabel('Filtered')
plt.grid(True)

plt.subplot(4, 1, 4)
plt.fill_between(np.arange(1,N), 0, lam)
plt.ylabel('Weights')
plt.grid(True)

plt.show()

