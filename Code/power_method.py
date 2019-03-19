from scipy.sparse.linalg import LinearOperator
from scipy.sparse import spmatrix
import numpy as np
import logging

def spectral_norm(A, tol=1e-8, max_iter=1000):
    """Computes the spectral norm of a linear operator A using power iteration.

    Parameters
    ===================
    - `A` (`numpy.ndarray`, `scipy.sparse.spmatrix`, or `scipy.sparse.linalg.LinearOperator`):
      the matrix for which we want to compute the spectral norm.

    Keyword parameters
    ====================
    - `tol` (float, default = `1e-8`): tolerance used to determine whether or not we
      should stop iterating. Once the estimates for the spectral norm are within distance
      `tol` of one another, we stop the power iterations and return.
    - `max_iter` (int, default = `1000`): maximum number of power iterations to do. If
      we reach this number of iterations then this function will return, but will display
      a warning that we reached the maximum number of iterations.
      - Power iteration can be extremely slow to converge, so you may need a large value
        of `max_iter` in order to find the true spectral norm.

    Return
    ====================
    - `sp_norm` (float): the estimated spectral norm of `A`.

    Code by Will Shand at the request of Stephen Becker, March 2019
    """
    if not any(issubclass(type(A),T) for T in [np.ndarray, spmatrix, LinearOperator]):
        raise ValueError("spectral_norm can only take arguments of type "
                "numpy.ndarray, scipy.sparse.spmatrix, or "
                "scipy.sparse.linalg.LinearOperator.")

    # Create an anonymous function matvec_op whose effect is equivalent to multiplying
    # the input by A'A.
    if issubclass(type(A), LinearOperator):
        matvec_op = lambda x: A.adjoint().matvec(A.matvec(x))
    else:
        matvec_op = lambda x: A.T.dot(A.dot(x))

    sp_norm = 0.
    sp_iter = np.random.normal(size = A.shape[-1])
    for ii in range(max_iter):
        Ax = matvec_op(sp_iter)
        new_sp_norm = np.linalg.norm(sp_iter)

        # Stopping condition when eigenvalue estimates get sufficiently close
        if abs(new_sp_norm - sp_norm) < tol:
            break
        else:
            sp_norm = new_sp_norm
            sp_iter = Ax / new_sp_norm

    if ii == max_iter-1:
        logging.warn(" spectral_norm ran for max_iter = %d iterations "
            "without converging. Returning..." % max_iter)

    return np.sqrt(sp_norm)

"""
TESTING
"""
if __name__ == "__main__":
    from scipy.sparse import random as sprandom

    # 1. Test on some random numpy arrays
    for ii in range(50):
        X = np.random.normal(size=(50,30))
        assert(abs(np.linalg.norm(X,2) - spectral_norm(X, max_iter=5000, tol=1e-8)) <= 1e-7)

    # 2. Test on some LinearOperator instances
    for ii in range(50):
        X1 = sprandom(50,50,density=0.2)
        X2 = np.random.normal(size=(10,50))

        # Dense representation of difference X1 - X2'X2
        A = X1 - X2.T.dot(X2)

        # LinearOperator representation of X1 - X2'X2
        mv  = lambda x: X1.dot(x) - X2.T.dot(X2.dot(x))
        rmv = lambda x: X1.T.dot(x) - X2.T.dot(X2.dot(x))
        L   = LinearOperator(X1.shape, matvec=mv, rmatvec=rmv)

        assert(abs(np.linalg.norm(A,2) - spectral_norm(L, max_iter=5000, tol=1e-8)) <= 1e-7)
