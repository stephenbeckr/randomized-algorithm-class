# Demos

Below are the demos from Spring 2019, which we will use and modify Fall 2021. Most are Matlab, some are Python or Julia.  TODO: create ipynb versions (with colab links) for the Python demos. (If a student wants to do this, please go ahead, and make a pull request)

The hyperlinks below are to ipynb (jupyter) notebooks using [nbviewer](https://nbviewer.jupyter.org) since github's default markdown interpreter doesn't always work (it often works if you refresh the page a few times, but not always)

- [Demo 1](https://nbviewer.jupyter.org/github/stephenbeckr/randomized-algorithm-class/blob/master/Demos/demo01_exactRankR.ipynb): a simple randomized SVD algorithm assuming exactly low-rank matrix, ([colab link for ipynb](https://colab.research.google.com/github/stephenbeckr/randomized-algorithm-class/blob/master/Demos/demo01_exactRankR.ipynb))
- [Demo 2](https://nbviewer.jupyter.org/github/stephenbeckr/randomized-algorithm-class/blob/master/Demos/demo02_sorts.ipynb): compare deterministic and randomized bubble and quicksorts, ([colab link for ipynb](https://colab.research.google.com/github/stephenbeckr/randomized-algorithm-class/blob/master/Demos/demo02_sorts.ipynb))
- Demo 3: compare different ways to compute the Frobenius norm of a matrix (in C)
- Demo 4: same as demo 3 but for sparse matrices (in Matlab). Comparse row vs column access
- Demo 5: speed/timing results for various Fast Johnson-Lindenstrauss Transforms
- Demo 6: statistical leverage scores applied to least-squares regressoin
- Demo 7: random matrix mulitplication via sub-sampling, following Drineas/Mahoney summer school notes
- Demo 8: high accuracy l2 regression via either iterative Hessian sketch or preconditioning (BLENDENPIK/LSRN)
- Demo 9: Randomized Kacmarz for solving consistent over-determined systems of equations
- Demo 10: l1 regression and p-stable distributions (for p=1,2, i.e., Cauchy and Gaussian)
- Demo 11: James-Stein estimator
- Demo 12: Basic noiseless Compressed Sensing demo
- Demo 13: Euclidean Distance Matrix (EDM) completion example, using nuclear norm minimization
- Demo 14: Monte Carlo integration and improvements (quasi-Monte Carlo, control variates), and comparison with quadrature
- Demo 15: Stochastic Gradient Descent (SGD) and improvements (SAGA, SVRG, Mini-batches, iterate averaging)
- Demo 16: Locality Sensitive Hashing (LSH): MinHash, SimHash, Euclidean distance hash
- Demo 17: LSH applied to k-Nearest-Neighbors (kNN)
- Demo 18: CountMin sketch to efficiently find frequencies of data names (data from Soc Security administration)
- Demo 19: AMS sketch vs Count sketch (median vs mean postprocessing)
- Demo 20: Core sets for Kmeans using Kmeans++ as coarse approximation

At some point I thought I had a distinction between "vignettes" and "demos" but I think that is gone now, and I've tried to rename them all to "demo"


### ipynb notebooks not rendering
Do you ever get the error message "Sorry, something went wrong. Reload?" when clicking on a `ipynb` file? If so, try refreshing the page a few times. If that doesn't resolve it soon (sometimes it does, sometimes it doesn't), then you can try either
1. go to <https://colab.research.google.com/> directly and "open" the file using the github interface
2. view the file by going to <https://nbviewer.jupyter.org/> and pasting in the URL
