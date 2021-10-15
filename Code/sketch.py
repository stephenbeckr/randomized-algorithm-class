"""
Implements some randomized linear sketches 
    (Gaussian, Haar, Count, FJLT with DCT, FJLT with Hadamard)
as well as some helper routines
    (Implicit2Explicit, TestAdjoints, TestSketch)

Part of APPM 5650 Randomized Algorithms

Taught/written by 
Stephen Becker, Oct 2021
stephen.becker@colorado.edu

It would be nice for someone who is competent with Python classes/inheritence 
to make this a proper subclass of LinearOperator and then add the methods,
 like Implicit2Explicit (as .todense or .toarray ) and TestAdjoints and
 TestSketch.

Major TODO: implement efficient adjoints for these operators (for the implicit ones,
 i.e., the FJLT ones, since the dense/sparse matrix based ones automatically
 have nice adjoints)

"""
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import LinearOperator, aslinearoperator, onenormest
from scipy.fft    import dct, idct

__all__ = ['fwht','Sketch','Gaussian','Haar','Count','FJLT','FJLT_Hadamard',
    'Implicit2Explicit','TestAdjoints','TestSketch']


def Sketch( fcnName, *args, **kwargs ):
    """
    Convenience function that wraps the actual sketching functions.
    Syntax of the form Sketch( sketchType, sz, rng )
        where sketchType is one of 'Gaussian','Haar','Count','FJLT','FJLT_Hadamard'
            (not case sensitive)
        and sz = (m,M) indicates the sketch S is of size m x M,
        and rng (optional) is of type np.random.default_rng

        Sketches are used to decrease dimensionality, so m <= M.  If m > M, 
        behavior of these functions is untested and may break.
    """

    sketchTypes = {'gaussian':Gaussian, 'jlt':Gaussian, 'haar':Haar,
        'count':Count, 'fjlt':FJLT, 'fjlt_dct':FJLT, 'fjlt_hadamard':FJLT_Hadamard }
    
    try:
        fcn=sketchTypes[fcnName.lower()]
    except KeyError:
        raise ValueError("Invalid sketch type: should be one of 'Gaussian','Haar','Count','FJLT','FJLT_Hadamard'") 
    return fcn( *args, **kwargs )

def Gaussian(sz, rng=np.random.default_rng() ):
    m, M    = sz
    return aslinearoperator( rng.standard_normal( size=(m,M) )/np.sqrt(m) )

def Haar(sz,  rng=np.random.default_rng() ):
    m, M    = sz
    # see also from scipy.stats import ortho_group  # Requires version 0.18 of scipy"
    A = rng.standard_normal( size=(M,m) )
    Q,R = np.linalg.qr(A, mode='reduced')
    sgn = np.sign( np.diag(R) )  # need to modify if A is complex...
    # we want A = (Q*sgn)*(sgn*R) in order to get a unique decomposition
    Q   = Q*sgn   # element-wise multiplication
    return aslinearoperator( Q.T*np.sqrt(M/m) )

def Count(sz,  rng=np.random.default_rng() ):
    # using ideas from https://github.com/scipy/scipy/blob/v1.7.1/scipy/linalg/_sketches.py
    #   written by Jordi Montes, 2017
    m, M    = sz
    rows = rng.integers(m, size=M, dtype = np.int64)
    #rows = rng_integers(rng, 0, m, M) # sklearn utility
    cols = np.arange(M+1)
    #signs = np.sign( rng.standard_normal(size=M) ).astype( np.int64 )
    signs = rng.choice([1, -1], M)
    S = scipy.sparse.csc_matrix((signs, rows, cols),shape=(m, M))
    return aslinearoperator(S)

def FJLT(sz,  rng=np.random.default_rng() ):
    m, M    = sz
    d       = np.sign( rng.standard_normal(size=M) ).astype( np.int64 ) # or rng.choice([1, -1], M)
    ind     = rng.choice( M, size=m, replace=False, shuffle=False)
    # IMPORTANT: make sure axis=0
    fjltMatMat = lambda X : np.sqrt(M/m)*_subsample( dct( _elementwiseMultiply(d,X), norm='ortho',type=3, axis=0) , ind)
    S       = LinearOperator( (m,M), matvec = fjltMatMat, matmat = fjltMatMat )
    return S

def FJLT_Hadamard(sz,  rng=np.random.default_rng() ):

    m, M    = sz
    M2       = _next_power_of_two(M)

    d       = np.sign( rng.standard_normal(size=M) ).astype( np.int64 )
    ind     = rng.choice( M2, size=m, replace=False, shuffle=False)
    fjltMatMat = lambda X : np.sqrt(1/m)*_subsample( fwht( _elementwiseMultiply(d,X)) , ind)
    S       = LinearOperator( (m,M), matvec = fjltMatMat, matmat = fjltMatMat )
    return S







def Implicit2Explicit( linOp, makeSparse = False, sparseFormat = 'csc' ):
    """ returns the explicit matrix representation of a linear operator 
        If makeSparse is True, then returns a sparse format
    """
    if not isinstance(linOp, LinearOperator) :
        raise ValueError('input must be a LinearOperator')
    
    if hasattr(linOp,'A'):
        A = linOp.A  # simple!
    else:
        m,n = linOp.shape
        A   = linOp@np.eye(n) 
    if makeSparse:
        if sparseFormat.lower() == 'csc':
            return scipy.sparse.csc_matrix(A)
        elif sparseFormat.lower() == 'csr':
            return scipy.sparse.csr_matrix(A)
        else:
            raise ValueError('Only sparse formats "csc" or "csr" are handled for now')
    else:
        return A
    # np.hstack([self.matvec(col.reshape(-1,1)) for col in X.T])

def TestAdjoints( A, At, method=None, rng = np.random.default_rng(), nReps = 10, 
        printMatrices = True, tol = 1e-10 ):
    if A.ndim != 2 or At.ndim != 2:
        raise ValueError("A and At must both be matrices ")
    shapeA  = A.shape
    shapeAt = At.shape
    if shapeAt[::-1] != shapeA:
        print("A and At do not have the same shape:", shapeA, shapeAt )
    m,n = shapeA

    if method is None:
        if m*n <= 100**2:
            method = 'explicit'
        else:
            method = 'implicit'
    
    if method.lower() == 'implicit':
        err = 0.
        for rep in range(nReps):
            x = rng.standard_normal( size=(n,1) )
            y = rng.standard_normal( size=(m,1) )
            err = np.max( [err,np.abs( y.T@(A*x) - (At@y).T@x )])

        normA = onenormest(A)  # spectral norm would be more classic, but this works
        normAt= onenormest(At)
        print(f"After {nReps:%d} trials, max error in inner product is {err:.2e}")
        print(f"  and ||A||_1 is {normA:.4e} while ||At||_1 is {normAt:.4e}")
    else:
        AA    = Implicit2Explicit(A)
        AAt   = Implicit2Explicit(At)
        err   = norm(AA-AAt.T)
        normA = norm(AA)
        normAt= norm(AAt)
        print(f"||A - At^T||_F is {err:.2e}")
        print(f"  and ||A||_F is {normA:.4e} while ||At||_F is {normAt:.4e}")
        if printMatrices:
            print("A is:")
            np.set_printoptions(precision = 2)
            print(AA)
            print("... and At is:")
            print(AAt)

    if err < tol:
        print("Passed check! These are likely adjoints")
        return True
    else:
        print("Failed check. Perhaps there is a bug")
        return False
    

def TestSketch( sz, style, nReps = 1000, printEvery = 100, rng = np.random.default_rng()  ):
    nReps = int(nReps)
    printEvery = int(printEvery)
    m,M = sz
    m   = int(m)
    M   = int(M)

    sumS = np.zeros( (M,M) )

    print('')
    print('Testing sketch of type', style)
    errList = []
    for rep in range(nReps):
        #S   = Sketch( (m,M), style, rng)
        S   = Sketch( style, (m,M), rng)
        if hasattr(S,'A'):
            A = S.A
        else:
            A = Implicit2Explicit( S )
            setattr(S,'A',A) # optional
            #raise ValueError('Need to implement implicit2explicit')
        sumS += A.T @ A
        errList.append( np.linalg.norm( sumS/(rep+1) - np.eye(M) )/M )

        if rep==0 or ( (rep +1) % printEvery == 0 ):
            print(f'Iter {rep+1:5d}, error is {errList[-1]:.2e}')

    sumS /= nReps
    print('The first 5 x 5 block of the sample mean is')
    with np.printoptions(precision = 2):
        print( sumS[:5,:5] )
    mn = np.mean(np.diag(sumS))
    print(f'and the avg diagonal entry (should be 1) is {mn:.7f}')






def fwht( x ):
    """ applies the Hadamard transform to x. If x has more than one column,
    the transform is applied to each column. 
    Leading dimension must be a power of 2; if not, you should zero-pad. This function
    doesn't do zero-padding for you.
    This code is not necessarily implemented well in terms of data copies and reshapes
    (and was written based off a Matlab code that assumed column-major format)

    To test that this code is correct, compare with
        import scipy.linalg
        H = scipy.linalg.hadamard(n, dtype = np.float)
    Stephen Becker, 9/23/2021
    """
    sz = x.shape
    m = sz[0]
    if not _is_power_of_two(m):
        #raise ValueError("Leading dimension of input must be a power of 2 for the Hadamard transform")
        df = _next_power_of_two(m) - m 
        if x.ndim == 1:
            y = np.pad( x, (0,df) )
        else:
            y = np.pad( x, ((0,df),(0,0)) )
        m2= _next_power_of_two(m)
    else:
        m2 = m
        y = x.copy()
    logm = int(np.log2(m2))

    # Assumes x is a numpy array, should probably check that (TODO)

    #for logk in range( int(np.log2(m)) ):
    for logk in range( logm ):
        k  = 2**logk      # 1, 2, 4, ..., m/2

        # use .resize to make permanent, but not recommended.
        yy = y.reshape( (2*k,-1, *sz[1:]), order='F')  # Need to "unpack" the tuple
        # in the above, order='F' (not the default order='C') is *very* important
        tmp = yy[:k, ...].copy()  # see https://stackoverflow.com/a/12116854, "Ellipsis"
        yy[:k, ...] += yy[k:, ...]
        yy[k:, ...]  = tmp - yy[k:, ...]

        y = yy.reshape( (m2,*sz[1:]), order='F') 
    
    return y

  
# ============= Helpers ===============
def _elementwiseMultiply( d, X):
    """ like d*X aka np.multiply(d,X)
    except it handles the case when d is size (n,) and X is size (n,1)
    since then naively doing d*X does an outer product since numpy doesn't
    consider (n,) and (n,1) to be the same...  but we also want to allow
    for the case when X is size (n,)
    """
    if d.ndim == X.ndim:
        # Great
        y = d*X
    elif d.ndim == 1:
        y = d.reshape(-1,1) * X
    else:
        y = d * X.reshape(-1,1)
    
    return y

def _subsample( X, ind):
    """ like X[ind,:] but works in case X has size (n,) """
    if X.ndim == 1:
        y = X[ind]
    elif X.ndim == 2:
        y = X[ind,:]
    else:
        raise ValueError("Expected 1D or 2D array")
    return y

def _is_power_of_two( n ):
    """ bit manipulations, suggested via https://stackoverflow.com/a/57025941 """
    n = int(n)
    return (n & (n-1) == 0) and n != 0

def _next_power_of_two(n):
    """ from https://stackoverflow.com/a/34625865 """
    n = int(n)
    return 1<<(n-1).bit_length()

if __name__ == '__main__':

    sketchList = ['Gaussian','Haar','Count','FJLT','FJLT_Hadamard']

    for sketch in sketchList:
        TestSketch( (10,20), sketch, nReps = 1000, printEvery=250)