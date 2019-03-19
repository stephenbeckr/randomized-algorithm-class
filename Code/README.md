# Code useful for the randomized algorithm class

## Misc utilities

- [my_normest](my_normest.m) for estimating the spectral norm in Matlab (variant of Matlab's normest that allows for function handles)
- [power_method](power_method.py) for estimating the spectral norm in Python
- [pdist2_faster](pdist2_faster.m) for calculating all pairwise distances, similar to Matlabs pdist2 but faster often (this version uses simple matrix multiplies)

## Hadamard Transform
[Hadamard Transform Code](https://github.com/jeffeverett/hadamard-transform) written by Stephen Becker and Jeff Everett, much faster than Matlab's code. You can also use this [Hadamard_teaching_code.m](Hadamard_teaching_code.m) simple .m file, which is actually faster than Matlab's code sometimes, and it's much simpler, so you can get a better idea of how the fast Hadamard transform works

## CountSketch

For a python version, see [scipy.linalg.clarkson_woodruff_transform](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.clarkson_woodruff_transform.html#scipy.linalg.clarkson_woodruff_transform). I haven't used it myself so not sure how efficient it is.

[CountSketch.c](countSketch.c) is an implementation of CountSketch in Matlab's mex interface. You'll need to compile it. See the source code for some fancy compilation options, but a basic version is:
```
>> mex countSketch.c
```
and a complicated faster version is:
```
>> mex countSketch.c -output countSketch_BLAS -DUSE_BLAS -lmwblas CFLAGS="\$CFLAGS -O3 -malign-double -march=native"
``` 
and run it like
```
>> mSmall = 10; mBig = 100; n = 7;
>> A = randn(mBig,n);
>> indx_map    = int64(randi(mSmall,mBig,1));
>> D = spdiags( sign(randn(mBig,1)), 0, mBig, mBig );
>> P = countSketch( D*A,  indx_map, mSmall, false );
>> P2 = countSketch( A'*D,  indx_map, mSmall, true )'; % alternative way
>> P2 = countSketch_BLAS( A'*D,  indx_map, mSmall, true )'; % if you did the fancy compile
```

[CountSketch_sparse.c](countSketch_sparse.c) is similar but works with sparse matrices. It doesn't have an option to use BLAS, so compilation is easy:
```
>> mex countSketch_sparse.c
```
This always assumes the transpose version. So, use it like:
```
>> A = sparse(A);
>> P2 = countSketch_sparse( A'*D,  indx_map, mSmall )';
```

