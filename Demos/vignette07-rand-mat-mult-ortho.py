# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# Vignette #XX
# APPM 4720/5720 Randomized Algorithms, Spring 2019
#
# Demo to show how to approximate a matrix-matrix product using randomization based on
# the approach of P. Drineas and M. Mahoney. See section four of "Lectures on Randomized
# Linear Algebra" by Drineas and Mahoney for additional details and analysis.
#
# This script focuses on the case where U \in \reals^{n \times d} with n >> d is an
# orthogonal matrix (specifically its columns are orthogonal). We choose c so that the
# deviation between our randomized matrix and the one we're approximating is bounded with
# high probability. The details for can be found in the Drineas and Mahoney paper.
#
# Algorithm randomized matrix multiplication.
# given U \in \reals^{n \times d} with n >> d, an integer c (1 \le c \le n),
#   and a probability distribution p of length n.
# repeat for k = 1,\dots, c:
#       1. pick i_k \in {1,\dots,n} with P(i = k) = p_k iid with replacement.
#       2. set R(k,:) = 1/sqrt(c * p_{i_k}) * U(i_k,:).
# return R^T R = sum_{k=1}^c 1/(c * p_{i_k}) * U(:,i_k) U(i_k,:).
# -------------------------------------------------------------------------------------- #
import numpy as np                      # import numpy package
import time as time                     # import time package
import matplotlib.pyplot as plt         # import matplotlib package
np.set_printoptions(precision = 2)      # display only four digits
# -------------------------------------------------------------------------------------- #
def vignette_rand_mat_mult_ortho(n_sims = 100, n = 1000, d = 6):
    # n_sims = 100; n = 1000; d = 6         # choose d ~ log n
    np.random.seed(seed = 2)                # set seed for reproducibility
    # ensure parameters are integers
    n_sims = np.int(n_sims);  n = np.int(n);  d = np.int(d);
    beta = 1; epsilon = 0.9; delta = 0.1;
    c_big_n = np.int(np.ceil(96*d/(beta*epsilon**2) *
                np.log(96*d/(beta*epsilon**2*np.sqrt(delta)))))
    c_sm_n = np.int(np.ceil(10*d**2 / (beta*epsilon**2)))
    c = np.amin([c_big_n, c_sm_n])
    # storage for simulations
    fro_norms = np.zeros(n_sims)
    # generate "data" matrices
    U = np.random.normal(scale = 1 / np.sqrt(n), size = (n,d))  # isotropic columns
    U, R_qr = np.linalg.qr(U)
    # ell_2 norms: columns of U and rows of U
    col_norm_U = np.linalg.norm(U, axis = 0); row_norm_U = np.linalg.norm(U, axis = 1)
    # probabilities for sampling
    probs = beta * row_norm_U**2 / d    # nearly optimal probs
    probs = probs / sum(probs)          # undo roundoff---ensure they sum to 1
    # compute theoretical upper bounds on expected squared Frobenius norms
    upper_bound = d**2 / (c*beta)          # optimal probabilities
    # simulation
    for t in range(n_sims):
        # initialize / re-zero matrices
        RTR = np.zeros((d,d))
        for k in range(c):
            # step 1
            i = np.random.choice(a = np.arange(n), replace = True, p = probs)
            # calculate rescaling coefficients
            rescale = 1 / np.sqrt(c*probs[i])
            # step 2
            R = rescale * U[i,:]
            # compute random matrix product via outerproduct
            RTR = np.outer(R, R) + RTR
        # calculate Frobenius norms of actual vs. randomized
        fro_norms[t] = np.linalg.norm(np.eye(d) - RTR, ord = 'fro')
    # print comparisons averaged across the number of simulations
    print('||U^T U - R^T R||_F^2 =', np.round(np.mean(fro_norms**2), 4),
          'vs', np.round(upper_bound, 4))

vignette_rand_mat_mult_ortho()

# Question: If you see that the estimated squared Frobenius norm for the naive
#           probabilities is smaller than that for the optimal probabilities, what should
#           you do?
