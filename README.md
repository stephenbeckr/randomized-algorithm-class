# Randomized Algorithm Class, APPM/STAT 5650, Fall 2021
Randomized algorithm class at CU Boulder, Fall 2021, [Professor Becker](http://amath.colorado.edu/faculty/becker/)

NOTE: Fall 2021, this is APPM/STAT 5650.  These course materials were created Spring 2019 when the same class was taught (also by Stephen Becker) under the special topics designation APPM 4720/5720; we are slowly updating the materials for Fall 2021 (there are no major changes -- the spirit of the class is the same as it was in 2019).  An earlier version of the special topics course was taught by Prof. Gunnar Martinsson.

Course information can be found at the public [course website](https://www.colorado.edu/amath/appm-47205720-open-topicsrandomized-algorithms-spring-2019)

The actual topics we covered, and links to references, are on this [google sheet](https://docs.google.com/spreadsheets/d/1z2yT99o8nCiotU0OZbrmmk0kAjff5iUDhKo3fpRVORA/edit?usp=sharing).  See below for a high-level list of what we covered. There was no single textbook for the class (and no standard set of notes).

This git repo contains things like code demos used in class. Most of the code is in Matlab; if any students want to translate demos to other languages and then push them, just make a pull request
- [Demos](Demos/)
- [Homeworks](Homeworks/) (homework solutions and code are on the private Canvas website)
- [Class policies, etc., for Spring 2019](APPM4720_5720_Spr2019_Syllabus.pdf)  (this will be updated for 2021)

# Student projections
Here is a [PDF showing a brief summary of the end-of-semester student projects from Spring '19](SlideshowAllPresentations_4720Spr19_Randomized.pdf).  If you're interested in the Shazam-like song matching algorithm (using locality sensitive hashing), their code is at [Padraig's github site](https://github.com/Lysandr/minHash_Shazam)

# Syllabus
- Introduction
  - Exact rank-r randomized SVD; randomized sorting
- How a computer works
  - Vector architecture, BLAS, LAPACK, blocking, memory hierarchy, 1-pass methods, sparse formats, HDF5
  - Reference: [Viktor Eijkhout's book](http://pages.tacc.utexas.edu/~eijkhout/istc/istc.html)
- Linear algebra
  - Norms and related inequalities, linear operators and operator norms, inner product on matrices, SVD, spectral and Frobenius norms, unitary invariance, Eckardt-Young, QR decomp.
- Basic probability
  - Conditional probability, expectation, variance, moments, linearity, law of total probability/variance/expectation, Bayes; Boole/union bound, Markov, Chebyshev, Jensen
- Tail bounds
  - Hoeffding/Chernoff/Bernstein
  - Reference: [Vershynin's 2018 book](https://www.math.uci.edu/~rvershyn/papers/HDP-book/HDP-book.html)
- Sketches, Johnson-Lindenstrauss Lemma
  - Reference for our version of proof: [Kakade and Shakhnarovich's notes](http://ttic.uchicago.edu/~gregory/courses/LargeScaleLearning/)
- JL to subspace embeddings via epsilon nets
  - Reference: [Woodruff 2014, Thm 2.3 and nearby](https://arxiv.org/abs/1411.4357)
- Fast Johnson-Lindenstrauss a la Ailon and Chazelle
  - References for detailed topics: [Becker and Pourkamali](http://arxiv.org/abs/1511.00152), [BLENDENPIK paper](https://dl.acm.org/citation.cfm?id=1958633), [Error correcting codes paper](http://proceedings.mlr.press/v37/ubaru15.html), [James' implementation code](https://github.com/jamesfolberth/fast_methods_big_data_project)
- CountSketch, very sparse sketches
  - References: [Clarkson and Woodruff for Countsketch](https://doi.org/10.1145/3019134 2017), [Li, Hastie and Church '06 for very sparse random projections](https://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf)
- Entry-wise sampling
  - References: [Achlioptas, Karnin and Liberty '13 "Simple is Best"](https://pdfs.semanticscholar.org/aa64/b8fb3382e42f90ee93a1dd0c78f13833963f.pdf), and also [Arora, Hazan and Kale '06](https://scholar.google.com/scholar?cluster=15000094124799237010&hl=en&as_sdt=0,6), [Achlioptas,  Karnin and Liberty '13](https://scholar.google.com/scholar?cluster=5902035377314465212&hl=en&as_sdt=0,6)
- Leverage scores, approximate matrix multiplication
  - References: [Mahoney's 2011 monograph](https://scholar.google.com/scholar?cluster=17892460197857789799&hl=en&as_sdt=0,6) and [Mahoney and Drineas 2017 summer school notes](http://arxiv.org/abs/1712.08880)
- Least-squares
  - deterministic backgroun: normal eq'ns, ill-conditioning, 1-pass classical algorithms, Tikhonov, Moore-Penrose pseudo-inverse
  - via sketching, Iterative Hessian Sketch, Blendenpik/LSRN, and randomized Kaczmarz
  - Reference for IHS: [Pilanci and Wainwright](http://www.jmlr.org/papers/volume17/14-460/14-460.pdf)
  - Reference for Blendenpik and LSRN: [Blendenpik paper](https://epubs.siam.org/doi/abs/10.1137/090767911), [LSRN paper](https://epubs.siam.org/doi/abs/10.1137/120866580)
  - Reference for randomized Kaczmarz: [Strohmer and Vershynin '08](http://www.springerlink.com/index/10.1007/s00041-008-9030-4), and extensions for importance sampling and SGD: [Eldar and Needell 2011](https://arxiv.org/abs/1008.4397), [Needell and Tropp 2012](https://arxiv.org/abs/1208.3805), and [Needell, Srebro and Ward 2016](https://arxiv.org/abs/1310.5715)
- l1 regression
  - Cauchy sketch, p-stable random variables
  - Reference: [Woodruff's monograph, ch 3](https://arxiv.org/abs/1411.4357), and [Clarkson et al.'s fast Cauchy transform](http://epubs.siam.org/doi/10.1137/140963698)
- randomized SVD
  - algorithm, one-pass variants, proofs
  - Reference: [Halko, Martinsson and Tropp 2011](https://epubs.siam.org/doi/10.1137/090771806), more practical versions [Tropp, Yurtsever, Udell and Cevher, SIMAX 2017](https://epubs.siam.org/doi/abs/10.1137/17M1111590)
  - Reference for one-pass version: [Yu et al, IJCAI 2017](https://arxiv.org/abs/1704.07669)
- Compressed sensing
  - Overview and proofs, quick slightly sub-optimal proofs of RIP via Johnson-Lindenstrauss
  - Reference: our proofs followed [Rauhut's 2011 monograph](http://www.mathc.rwth-aachen.de/~rauhut/files/LinzRauhut.pdf)
- Matrix Completion/Euclidean Distance Completion
- Monte Carlo
  - Background, background on quadrature
  - Improvements to Monte Carlo (quasi-Monte Carlo and control variates)
  - Reference for quasi-MC: [Dick, Kuo, Sloan; Acta Numerica 2013](https://doi.org/10.1017/S0962492913000044)
- Stochastic Gradient Descent (SGD)
  - Improvements (SAGA, SVRG)
  - Reference: [Bottou, Curtis and Nocedal, Siam Review '18](http://arxiv.org/abs/1606.04838)
- Locality Sensitive Hashing (LSH)
  - Background on hash functions (hash tables, cryptographic hashes)
  - k-NN, MinHash, SimHash, Euclidean distance
  - Reference for LSH: [ch 3 of Rajaraman and Ullman 2010](http://infolab.stanford.edu/~ullman/mmds/ch3a.pdf)
  - Reference for hashing: A classic reference is vol 3 of D. Knuth's 1968 "The Art of Computer Programming", but wikipedia is also convenient
- CountMin sketch and friends
  - min/median non-linear post-processing; AMS sketch; versions of Count sketch
  - Reference: [Cormode's 2011 review](http://www.cs.umass.edu/~mcgregor/711S12/sketches1.pdf) and [Cormode's 2013 Simon's talk](http://dimacs.rutgers.edu/~graham/pubs/html/TalkSimons13.html)
- Coresets
  - For k-means clustering via k-means++
  - Reference: [Bachem, Lucic and Krause '17](https://arxiv.org/abs/1703.06476)

Guest lectures:
- Richard Border about [Stochastic Lanczos quadrature, Border and Becker '19](https://www.biorxiv.org/content/10.1101/607168v1)
- Osman Malik about tensor sketches and interpolative decomposition: [Malik and Becker 2018 NeurIPS](https://papers.nips.cc/paper/8213-low-rank-tucker-decomposition-of-large-tensors-using-tensorsketch) and [Malik and Becker 2019 ID](https://arxiv.org/abs/1901.10559)

# Resources (since there is no main textbook)
Here is a list of monographs that are helpful:

- Gunnar Martinsson previously taught a version of this course, and still has information on the course website [Fast algorithms for big data from spring 2016](https://amath.colorado.edu/faculty/martinss/Teaching/APPM5720_2016s/index.html).
- Gunnar Martinsson also has a 34 page monograph [Randomized methods for matrix computation (started 2016, updated 2018)](http://arxiv.org/abs/1607.01649)
- Roman Vershynin published a Nov 2018 300 page book [High-Dimensional Probability: An Introduction with Applications in Data Science](https://www.cambridge.org/gb/academic/subjects/statistics-probability/probability-theory-and-stochastic-processes/high-dimensional-probability-introduction-applications-data-science?format=HB#8FyVeUCeT4PcEjMp.97) published by Cambridge University Press
  - Like many books from Cambridge University Press, they are allowing Vershynin to host a [PDF copy on his website](https://www.math.uci.edu/~rvershyn/papers/HDP-book/HDP-book.html).
  - A good complement is [Joel Tropp's 2011 "An Introduction to Matrix Concentration Inequalities"](https://arxiv.org/abs/1501.01571)
- David P. Woodruff has a 157 page monograph [Sketching as a Tool for Numerical Linear Algebra](http://dx.doi.org/10.1561/0400000060) from 2014 as part of the Foundations and Trends in Theoretical Computer Science series.
  - You can access the text via the [arXiv version](https://arxiv.org/abs/1411.4357)
- Michael Mahoney has a 100 page monograph [Randomized algorithms for matrices and data](http://dx.doi.org/10.1561/2200000035) from 2011 as part of the Foundations and Trends® in Machine Learning series
  - You can access the text via the [arXiv version](https://arxiv.org/abs/1104.5557)
- [Lectures on Randomized Numerical Linear Algebra, by Michael Mahoney and Petros Drineas, 2017](http://arxiv.org/abs/1712.08880), "This chapter is based on lectures on Randomized Numerical Linear Algebra from the 2016 Park City Mathematics Institute summer school on The Mathematics of Data."
- Alex Gittens teaches [CSCI 6220 Randomized Algorithms](http://www.cs.rpi.edu/~gittea/teaching/fall2018/csci6220-and-4030.html) at RPI and covers good material

## Textbooks for related topics
For numerical linear algebra, see [Applied Numerical Linear Algebra by James Demmel (SIAM 1997)](https://epubs.siam.org/doi/book/10.1137/1.9781611971446?mobileUi=0). You have free access to a PDF version of each chapter if you are on CU campus
