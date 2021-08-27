#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "blas.h"

/* Stephen Becker, Jan 15 2018
 * Make sure to allocate the array A on the heap,
 * not on the stack!
 *  e.g.,  double A[1000*1000]; is not a good idea
 *
 * Usage:
 *   gcc demo03_FrobeniusNorm.c -O3
 *   ./a.out n      # uses a n x n matrix
 *   ./a.out n 1    # if a 2nd argument is positive, use row-based version
 *   gcc demo03_FrobeniusNorm.c -O3 -I/Applications/MATLAB_R2017b.app/extern/include -lblas
 *   ./a.out n 1 1   # if a 3rd argument is positve, use BLAS for rows/columns
 *      ... if this 3rd argument is negative, use BLAS and vectorize
 *
 * Results for 10,000 x 10,000:
 * NO BLAS:
 *  2.8 s   Looping over rows, inner loop over columns
 *  .75 s   Looping over columns, inner loop over rows
 * WITH BLAS:
 *  2.5 s   Looping over rows, inner loop over columns
 *  .52     Looping over columns, inner loop over rows
 *  */

/* See dnrm2( ptrdiff_t *N, double *X, ptrdiff_t *INCX)  */

int main(int argc, char *argv[]) {

    ptrdiff_t m, n, length;
    ptrdiff_t INCX;
    int     i, j; /* counters */
    int     ROW_BASED, USE_BLAS, VECTORIZE = 0; /* boolean flags */
    double  *A;
    double  s=0.; /* sum */
    double  t = 0.;


    m   = 10;
    if (argc > 1)
        m = atoi( argv[1] );
    n   = m;


    /* We're storing it in column-major format */
    /* The problem is, this code can be slow itself,
     *  so for speed runs, set it all to zero and hope
     *  compiler doesn't try to be too clever */
    /*
    A   = malloc( m * n * sizeof( double ) );
    for (i=0; i<m; i++ )
        for (j=0; j<n; j++ )
            A[i + m*j] = i + m*j;
    */

    A   = calloc( m * n , sizeof( double ) ); 

    printf("A[1] is %f and A[end] is %f\n", A[1], A[ m*n - 1] );


    if ( (argc > 2 ) && (atoi(argv[2])>0) )
        ROW_BASED = 1;
    else
        ROW_BASED = 0;
    if ( (argc > 3 ) && (atoi(argv[3])>0) )
        USE_BLAS = 1;
    else if ( (argc > 3 ) && (atoi(argv[3])<0) ) {
        USE_BLAS = 1;
        VECTORIZE = 1;
    }    
    else
        USE_BLAS = 0;
    if (USE_BLAS)
        printf("Using BLAS\n");
    if (VECTORIZE){
        printf("Vectorizing (this is the *proper* way, no 'for' loops)\n");
        INCX = 1;
        length = m + n;
        s = dnrm2( &length, A, &INCX );

        free( A );
        printf("... sum of squared entries is %e\n", s );

        return 0;


    }

    if (ROW_BASED == 1 ){
        printf("Looping over columns, inner loop over rows\n");
        if (USE_BLAS) {
            INCX = 1;
            for (j=0; j<n; j++) {
                t   = dnrm2( &m, A + m*j, &INCX );
                s   += t*t;
            }
        } else {
            for (j=0; j<n; j++)
                for (i=0; i<m; i++ )
                    s += A[i+m*j]*A[i+m*j];
        }
    } else {
        printf("Looping over rows, inner loop over columns\n");
        if (USE_BLAS) {
            INCX = m;
            for (i=0; i<m; i++) {
                t   = dnrm2( &n, A+i, &INCX );
                s   += t*t;
            }

        } else {
            for (i=0; i<m; i++ )
                for (j=0; j<n; j++)
                    s += A[i+m*j]*A[i+m*j];
        }
    }
    free( A );
    printf("... sum of squared entries is %e\n", s );

    return 0;

}
