# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# Vignette #XX
# APPM 4720/5720 Randomized Algorithms, Spring 2019
#
# Demo to show how to approximate a matrix-matrix product using randomization based on
# the approach of P. Drineas and M. Mahoney. See section four of "Lectures on Randomized
# Linear Algebra" by Drineas and Mahoney for additional details and analysis.
#
# Algorithm randomized matrix multiplication.
# given A \in \reals^{m \times n}, B \in \reals^{n \times p}, an integer c
#    (1 \le c \le n), and a probability distribution p of length n.
# repeat for k = 1,\dots, c:
#       1. pick i \in {1,\dots,n} with P(i = k) = p_k iid with replacement.
#       2. set C(:,k) = 1/sqrt(c * p_i) * A(:,i) and R(k,:) = 1/sqrt(c * p_i) * B(i,:).
# return C, R, and CR = sum_{k=1}^c 1/(c * p_i) * A(:,i) B(i,:).
# -------------------------------------------------------------------------------------- #
import numpy as np                      # import numpy package
import time as time                     # import time package
import matplotlib.pyplot as plt         # import matplotlib package
np.set_printoptions(precision = 2)      # display only four digits
# -------------------------------------------------------------------------------------- #
def vignette_rand_mat_mult(n_sims = 1000, m = 100, n = 20, p = 80):
    np.random.seed(seed = 2)                # set seed for reproducibility
    # ensure parameters are integers
    n_sims = np.int(n_sims);  m = np.int(m);  n = np.int(n);  p = np.int(p)
    c = min(max(np.int(np.round(0.5*n)), 1), n)
    # storage for simulations
    fro_norms = np.zeros(n_sims);   fro_norms_opt = np.zeros(n_sims)
    # generate "data" matrices
    A = np.random.normal(scale = 1 / np.sqrt(n), size = (m,n))  # isotropic rows
    B = np.random.normal(scale = 1 / np.sqrt(n), size = (n,p))  # isotropic columns
    AB = np.matmul(A, B);
    # ell_2 norms: columns of A and rows of B
    col_norm_A = np.linalg.norm(A, axis = 0); row_norm_B = np.linalg.norm(B, axis = 1)
    # probabilities for sampling
    probs = np.ones(n)/n                # naive probabilities---i.e., uniform distribution
    probs_opt = (col_norm_A * row_norm_B) / sum(col_norm_A * row_norm_B) # optimal probs
    probs_opt = probs_opt / sum(probs_opt)  # undo roundoff---ensure they sum to 1
    # compute theoretical upper bounds on expected squared Frobenius norms
    upper_bound = sum(col_norm_A**2 * row_norm_B**2 / (c*probs))   # naive probabilities
    upper_bound_opt = sum(col_norm_A * row_norm_B)**2 / c          # optimal probabilities
    # simulation
    for t in range(n_sims):
        # initialize / re-zero matrices
        C = np.zeros((m,c));     R = np.zeros((c,p))
        C_opt = np.zeros((m,c)); R_opt = np.zeros((c,p))
        for k in range(c):
            # step 1
            i = np.random.choice(a = np.arange(n), replace = True, p = probs)
            i_opt = np.random.choice(a = np.arange(n), replace = True, p = probs_opt)
            # calculate rescaling coefficients
            rescale = 1 / np.sqrt(c*probs[i]);
            rescale_opt = 1 / np.sqrt(c*probs_opt[i_opt])
            # step 2
            C[:,k] = rescale * A[:,i]; R[k,:] = rescale * B[i,:]
            C_opt[:,k] = rescale_opt * A[:,i_opt]; R_opt[k,:] = rescale_opt * B[i_opt,:]
        # compute random matrix product via outerproduct
        CR = np.matmul(C, R)
        CR_opt = np.matmul(C_opt, R_opt)
        # calculate Frobenius norms of actual vs. randomized
        fro_norms[t] = np.linalg.norm(AB - CR, ord = 'fro')**2
        fro_norms_opt[t] = np.linalg.norm(AB - CR_opt, ord = 'fro')**2
    # print comparisons averaged across the number of simulations
    print('||AB - CR||_F^2 simulation vs. upper bound:',
            np.round(np.mean(fro_norms),2), 'vs.',
            np.round(upper_bound,2), '(using naive sampling probabilities)')
    print('||AB - CR||_F^2 simulation vs. upper bound:',
            np.round(np.mean(fro_norms_opt),2), 'vs.',
            np.round(upper_bound_opt,2), '(using optimal sampling probabilities)')
    # return np.mean(fro_norms), upper_bound, np.mean(fro_norms_opt), upper_bound_opt

# small problem
vignette_rand_mat_mult(n_sims = 10000, m = 10, n = 4, p = 10)

# moderate problem
vignette_rand_mat_mult(n_sims = 1000, m = 1000, n = 20, p = 1000)

# big problem
vignette_rand_mat_mult(n_sims = 1000, m = 5000, n = 10, p = 5000)

# Question: If you see that the estimated squared Frobenius norm for the naive
#           probabilities is smaller than that for the optimal probabilities, what should
#           you do?
