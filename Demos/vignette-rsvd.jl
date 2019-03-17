#=
Algorithm RSVD.
given a matrix A \in \reals^{m \times n}, a target rank k and an oversampling
      parameter p (e.g., p = 10 is a good choice).
stage A. Find an approximate range.
    1. form an n \times (k+p) Gaussian random matrix G.
    2. form the sample matrix Y = AG.
    3. orthonormalize the columns of Y via a QR factorization.
stage B. Form a specific factorization.
    4. form the (k+p) \times n matrix B = Q'A.
    5. form the SVD of the (small) matrix B as B = \hat{U} D V'.
    6. form U = Q \hat{U}.
return matrices U, D, and V as an approximate rank (k+p) SVD of A.
=#
using Random, LinearAlgebra, Plots, Statistics
# ---------------------------------------------------------------------------- #
function vignette_rsvd()
rng = Random.seed!(2);              # set seed for reproducibility
n_sims = 10;                        # number of simulations
n_subs = 25;                        # number of subsamples
m = 2000;                           # rows of matrix
n = 20*ceil(log(m));                # columns of matrix
# n_sims = 50; n_subs = 25; m = 1000; n = 15*ceil(log(m)); # alt. run parameters
# n_sims = 1;  n_subs = 25; m = 5000; n = 15*ceil(log(m)); # alt. run parameters
k = ceil(log(m));                   # target rank
p = max(ceil(log(m)), 10);          # oversampling parameter
n = Int16(n); k = Int16(k); p = Int16(p);   # convert floats to integers
# ---------------------------------------------------------------------------- #
times = zeros(n_sims, 2);           # bookkeeping for run times
fro_mean = zeros(n_sims, n_subs);   # bookkeeping for mean Frobenius norm
op_mean = zeros(n_sims, n_subs);    # bookkeeping for mean operator norm
norm_bound = zeros(n_sims, 2);      # bookkeeping for theoretrical norm
fro_lo = zeros(n_sims, 1);          # bookkeeping for Frobenius norm bound
op_lo = zeros(n_sims, 1);           # bookkeeping for operator norm bound
# ---------------------------------------------------------------------------- #
for j in 1:n_sims
    # "data" matrix
    A = [2*ones(m,2)+rand(m,2) randn(m, k-2) 0.01*randn(m, n-k)]/sqrt(m);
    # NOTE: the matrix A has k columns that lead to "important" singular values.
    # The remaing n-k columns correspond to fast-decaying singular values.
    for i in 1:n_subs
        # stage A of randomized SVD
        G = randn(n, k+p);                              # Gaussian random matrix
        Y = A*G;                                        # sample matrix
        F = qr(Y); Q = Matrix(F.Q);                     # orthonormalize Y
        # stage B of randomized SVD
        B = Q'*A;                                       # form small matrix
        U_B, D_B, V_B = svd(B);                         # SVD of B
        U = Q*U_B;                                      # rank k matrix
        # bookkeeping and comparisons
        U_A, D_A, V_A = svd(A);                         # SVD of A
        fro_mean[j,i] = norm(A - Q*Q'*A, 2);
        op_mean[j,i] = opnorm(A - Q*Q'*A);
    end # end of inner simulation loop (i.e., simulation for a fixed matrix A)
    fro_lo[j,1] = sum( D_A[k+1:min(m,n)].^2 )^0.5;
    op_lo[j,1] =  D_A[k+1];
    norm_bound[j, 1] = (1 + k / (p-1))^(0.5) * fro_lo[j,1];     # Frobenius
    norm_bound[j, 2] = (1 + sqrt(k / (p-1))) * op_lo[j,1] +     # operator
                            exp(1) * (sqrt(k+p) / p) * fro_lo[j,1];
end # end of outer simulation loop
# ---------------------------------------------------------------------------- #
# plot showing singular value decay
p3 = plot(D_A, yscale = :log10, linecolor = :blue,
            marker = :circle, markercolor = :blue, label = "full",
            title = "sing. value decay (final sim.)");
p4 = plot(D_A[1:k+p], yscale = :log10, linecolor = :blue,
            marker = :circle, markercolor = :blue, label = "full",
            title = "sing. value comparison (final sim.)")
plot!(D_B, yscale = :log10, linecolor = :red, linestyle = :dash,
            marker = :x, markerstrokecolor = :red, label = "randomized");

# ---------------------------------------------------------------------------- #
# plots of average Frobenius and operator norms vs. theoretical bounds
p1 = plot(norm_bound[:,1],linecolor = :blue, marker = :circle,
            markercolor = :blue, label = "upper",
            title = "E-Y bounds: Frobenius");
plot!(mean(fro_mean, dims = 2),linecolor = :red, linestyle = :dash,
     marker = :x, markerstrokecolor = :red, label = "mean");
     plot!(fro_lo, linecolor = :blue, linestyle = :dot,
           marker = :star8, markercolor = :blue, label = "lower",
           legend = :bottomright);
p2 = plot(norm_bound[:,2],linecolor = :blue, marker = :circle,
            markercolor = :blue, label = "upper",
            title = "E-Y bounds: operator");
plot!(mean(op_mean, dims = 2),linecolor = :red, linestyle = :dash,
     marker = :x, markerstrokecolor = :red, label = "mean");
plot!(op_lo,linecolor = :blue, linestyle = :dot,
      marker = :star8, markercolor = :blue, label = "lower",
      legend = :bottomright);
# return summary plot as "output" of the function
plot(p3, p4, p1, p2, layout=(2,2))
end # end of function
# ---------------------------------------------------------------------------- #
vignette_rsvd()
