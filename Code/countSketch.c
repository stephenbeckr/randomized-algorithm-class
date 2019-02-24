/*
 * Standard Usage:
 *
 * indx_map    = int64(randi(mSmall,mBig,1));
 * D = spdiags( sign(randn(mBig,1)), 0, mBig, mBig );
 *
 * doTranspose = false;
 * P = countSketch( D*A,  indx_map, mSmall, doTranspose )
 *
 * or 
 *
 * doTranspose = true;
 * P  = countSketch( A'*D, indx_map, mSmall, doTranspose )'
 * (note the transposes; this version should be faster
 *  due to how Matlab stores matrices, and the fact that
 *  Matlab can transpose a matrix very efficiently)
 * 
 * A        is an mBig x N matrix
 * P        is an mSmall x N matrix
 * indx_map is an int64 vector of length mBig where each entry
 *      is an integer in [1,mSmall] (i.e., 1-based indexing, Matlab-style)
 *
 *
 * Note: the random sign flips are NOT done in this mex file,
 *  that is why you should multiply by the diagonal D matrix as suggested
 *  above. The code is written this way because Matlab is very efficient
 *  at the diagonal multiply, so writing that myself in C would lead to
 *  worse performance. Or maybe it would be fast when using the BLAS
 *  version but not by much, so I was too lazy to do it...
 *
 * Implements the "CountSketch"
 * as known in the data streaming literature
 * (e.g., 7] Moses Charikar, Kevin Chen, and Martin Farach-Colton. 
 * "Finding frequent items in data streams". Theor.
 *  Comput. Sci., 312(1):3–15, 2004 )
 *
 * In the compressed least squares literature, this was analyzed in
 *
 * "Low Rank Approximation and Regression in Input Sparsity Time"
 * Kenneth L. Clarkson, David P. Woodruff
 * http://arxiv.org/abs/1207.6365
 * STOC '13, co-winner of best paper award
 *
 * Computational complexity is nnz(A)
 *
 * */

/*  Compiling:
 *
 * Compile with just:  mex countSketch.c
 *
 * or you can try to get more performance with things like
 *
 * mex countSketch.c CFLAGS="\$CFLAGS -O3 -malign-double -march=native" 
 *
 * or 
 *
 * mex countSketch.c -output countSketch_BLAS -DUSE_BLAS -lmwblas CFLAGS="\$CFLAGS -O3 -malign-double -march=native" 
 *
 * (these are for GCC compilers, you'll need to change the flags slightly
 *  for MVCC or LLVM/Clang compilers, e.g., Clang doens't use -malign-double )
 * NOTE: pass the flag -DLONGLONG to use LONG LONG pointer to indx_map. This is 
 * necessary if LONG is of size 32 bit.
 *
 */

/* Notes:
 *
 * For efficiency, since Matlab uses column-major order,
 * the input should be At ( = A') and NOT A
 * Likewise, the output is Pt ( = P' = (Sketch(A))' )
 *
 * In theory, this can be applied to sparse matrices
 * It would only be efficient if they are stored in csr order
 * (Matlab uses csc), or if we have the transpose of a csc matrix
 * (i.e., do the exact same transpose trick we do for sparse
 *  matrices )
 *
 * For now, does NOT do sparse matrices
 * There is a separate code just for sparse matrices
 *
 * Warning: the code does not do error checks, so it can easily crash.
 * Make sure that the "indx_map" is of type int64
 *
 * Stephen Becker, srbecker@us.ibm.com, June 5 2014
 * The use of CountSketch was suggested by Haim Avron
 *
 * Updates by Stephen Becker, Feb 2019
 **/
#if defined(__GNUC__) && !(defined(__clang__)) 
#include <uchar.h>
#endif
#include <math.h>
#include "mex.h"

#ifdef USE_BLAS
#include "blas.h"
#endif


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *At,*Pt; 
    double *A, *P;
    double alpha;
#ifdef LONGLONG
    long long  *indx_map;
#else
    long   *indx_map;
#endif
    mwSize mBig,mSmall,n, i,j,k;
    int     DO_TRANSPOSE=1;
#ifdef USE_BLAS
    ptrdiff_t   size, stride, strideBig;
#endif
    
    /* Check for proper number of arguments */
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "MATLAB:countSketch:invalidNumInputs",
                "Four input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:countSketch:maxlhs",
                "Too many output arguments."); 
    } 
    if ( !(mxIsInt64(prhs[1])) )
        mexErrMsgTxt("2nd input must be of type int64");
#ifdef LONGLONG
    indx_map      = (long long *)mxGetData( prhs[1] );
#else
    indx_map      = (long *)mxGetData( prhs[1] );
#endif
    mSmall        = mxGetScalar( prhs[2] );
    DO_TRANSPOSE  = (int)mxGetScalar( prhs[3] );
    
    if (mxIsSparse(prhs[0]))
        mexErrMsgTxt("Cannot handle sparse 'A' matrix, try countSketch_sparse.c instead");
    
    if ( DO_TRANSPOSE == 1 ) {
        At  = mxGetPr(prhs[0] );
        n   = mxGetM( prhs[0] );
        mBig= mxGetN( prhs[0] );
        /* Create a matrix for the return argument */
        plhs[0] = mxCreateDoubleMatrix( (mwSize)n, (mwSize)mSmall, mxREAL);
        Pt      = mxGetPr( plhs[0] );
        P       = NULL; /* try to prevent bugs */
        A       = NULL;
#ifdef USE_BLAS
        size    = (ptrdiff_t)n;
        stride  = (ptrdiff_t)1;
#endif
        /* And the actual computation:
         * Copy columns of At to Pt */
        alpha = 1.;
        for ( i=0; i < mBig ; i++ ) {
            k   = indx_map[i]-1; /* 0-based */
            /* copy Pt(:,k) <-- At(:,i)
             * e.g. since height of Pt is N,
             * P + k*n <-- At + i*n  */
#ifdef USE_BLAS
            daxpy(&size,&alpha,At+i*n,&stride,Pt+k*n,&stride);
#else
            for ( j=0; j<n; j++ )
                Pt[k*n+j] += At[i*n+j];
#endif       
        }
        
    } else if ( DO_TRANSPOSE == 0 ) {
        A   = mxGetPr(prhs[0] );
        n   = mxGetN( prhs[0] );
        mBig= mxGetM( prhs[0] );
        /* Create a matrix for the return argument */
        plhs[0] = mxCreateDoubleMatrix( (mwSize)mSmall, (mwSize)n, mxREAL);
        P       = mxGetPr( plhs[0] );
        Pt      = NULL;
        At      = NULL;
#ifdef USE_BLAS
        size    = (ptrdiff_t)n;
        stride  = (ptrdiff_t)mSmall;
        strideBig  = (ptrdiff_t)mBig;
#endif
        
        /* And the actual computation:
         * Copy rows of A to P */
        alpha = 1.;
        for ( i=0; i < mBig ; i++ ) {
            k   = indx_map[i]-1; /* 0-based */
            /* copy P(k,:) <-- A(i,:) */
#ifdef USE_BLAS
            daxpy(&size,&alpha,A+i,&strideBig,P+k,&stride);
#else
            for ( j=0; j<n; j++ )
                P[k+j*mSmall] += A[i+j*mBig];
#endif
            
        }
        
    } else {
        mexErrMsgTxt("4th input must 0 (no transpose) or 1 (transpose)");
    }


    
    
    return;
}
