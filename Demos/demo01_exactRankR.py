# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# Vignette #1
# APPM 4720/5720 Randomized Algorithms, Spring 2019
# Stephen Becker (original MATLAB) & Jake Knigge (Python)
# This is not a practical algorithm, since it only works if the matrix is *exactly* low
# rank.
# -------------------------------------------------------------------------------------- #
import numpy as np                      # import NumPy package
import time as time                     # import time package
# -------------------------------------------------------------------------------------- #
np.set_printoptions(precision = 4)      # display only four digits
np.random.seed(seed = 2)                # set seed for reproducibility
n = np.int(4e3); m = n                  # dimension of problem
r = np.int(100)                         # rank of matrix
mu, sigma = 0, 1                        # mean and standard deviation
zz = np.random.normal(mu, sigma, n*r)   # generate (normal) random numbers
Z = zz.reshape(n,r)                     # reshape to matrix
A = np.matmul(Z, Z.T)                   # compute outerproduct matrix
# -------------------------------------------------------------------------------------- #
# Find its SVD with conventional methods
t = time.time()                         # time SVD calculation
U, S, Vh = np.linalg.svd(A); V = Vh.T   # compute SVD of A and transpose V
elapsed = time.time() - t
print('The full SVD took', round(elapsed, 4), 'seconds.')
# -------------------------------------------------------------------------------------- #
# Find its SVD with a randomized method
tt = time.time()
t = time.time(); Omega = np.random.normal(mu, sigma, (n, r)); print(round(time.time() - t, 4), 'seconds')
t = time.time(); Y     = np.matmul(A, Omega);                 print(round(time.time() - t, 4), 'seconds')
t = time.time(); Q, R  = np.linalg.qr(Y, mode='reduced');     print(round(time.time() - t, 4), 'seconds')
t = time.time(); QtA   = np.matmul(Q.T, A);                   print(round(time.time() - t, 4), 'seconds')
tm  = time.time() - tt
print('The approximate SVD took', round(tm, 4), 'seconds.')
# -------------------------------------------------------------------------------------- #
A_estimate = np.matmul(Q, QtA)
err = np.linalg.norm( np.ravel(A) - np.ravel(A_estimate) ) / np.linalg.norm( np.ravel(A) )
print('The error ||A-A_estimate||_F/||A||_F is ', '{:0.4e}'.format(err), '.', sep = '')
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
