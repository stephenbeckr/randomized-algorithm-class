# ---------------------------------------------------------------------------- #
#=
 Testing computational speed of various Fast JL transforms
  Feb 11 2019
Y = S*X, X is M x N and S is m x M

Code by: Stephen Becker
Julia modifications: Jake Knigge
=#
# ---------------------------------------------------------------------------- #
# load packages
using Random, LinearAlgebra, Plots, Statistics, Hadamard, FFTW, SparseArrays
# ---------------------------------------------------------------------------- #
function Hadamard_teaching_code(x)
# Hadamard_teaching_code(x)
#   applies the Hadamard transform to x.
#   If x has more than one column, the transform is applied
#   to each column.
# This code is not fast, but it shows you how to exploit
#   the structure of the transform.
# Note: this code does not do any sub-sampling.
# Reference: https://en.wikipedia.org/wiki/Fast_Walsh–Hadamard_transform
  m = size(x,1); n = size(x,2);
  if nextpow(2,m) != m
      print("ERROR: The number of rows of x must be a power of 2.")
      return
  end
  y = copy(x)
    for bit in 1:log2(m)
        k   = Int64(2^bit);      # e.g., 2, 4, ..., m
        k2  = Int64(2^(bit-1));  # e.g., 1, 2, ..., m/2
        y   = reshape(y, k, :, n);
        tmp = y[1:k2,:,:];
        y[1:k2,:,:]       = y[1:k2,:,:] + y[k2+1:k,:,:];
        y[(k2+1):k,:,:]   = tmp         - y[k2+1:k,:,:];
        y   = reshape(y, m, n)
    end # loop
    return y
end # function
# ---------------------------------------------------------------------------- #
function slowCountSketch( DX, targetRows )
# slow version of count sketch
    m   = length( targetRows );
    Y   = zeros(m, size(DX,2) );
    for j in 1:size(DX,1)
        i  = targetRows[j];
        Y[i,:] = Y[i,:] + DX[j,:];
    end # loop
end # function
# ---------------------------------------------------------------------------- #
# test to generate an error
x = randn(65,5);
y = Hadamard_teaching_code(x)

# test to generate an error
x = randn(64,5);
y = Hadamard_teaching_code(x);
y_alt = 64*fwht_natural(x, 1);
norm(y - y_alt)
# Check Implementations for correctness
X = randn(2^13,100); M, N = size(X); m = M/4;

# == Hadamard Code ==
# Check for normalization and ordering
Y1 = Hadamard_teaching_code(Matrix{Float64}(I, 4, 4))
Y2 = 4*fwht_natural(Matrix{Float64}(I, 4, 4), 1)
norm(Y1 - Y2)
# Compare times---run more than once b/c of Julia's "just-in-time" compiler!
@time Y1 = Hadamard_teaching_code(X);
@time Y2 = M*fwht_natural(X, 1);
norm(Y1 - Y2)
# ---------------------------------------------------------------------------- #
# Test speed
N = 100;
M_list  = 2 .^(11:14)# (^).(2,10:13); # defined using "broadcast" operation
nTrials     = 10; # get some averages
nAlgos      = 6; # change to 7 if/when CountSketch is included
Times       = zeros(nAlgos,length(M_list),nTrials);
Times_setup = Times;
ALGO_NAMES = ["Gaussian","FJLT - DCT","FJLT - Hadamard", #"Count",
              "Very Very sparse","Very sparse","Sparse"];
rng = Random.seed!(9);
# ---------------------------------------------------------------------------- #
for Mi in 1:length( M_list )
  println("Mi is ", Mi, " of ", length(M_list), ".");
  for trial = 1:nTrials
    M   = M_list[Mi];
    m   = round(M/4);
    X   = randn(M,N);

    ALGO = 1; # Gaussian sketch
    tic = time();
    S = randn(Int64(m),M);
    Times_setup[ALGO,Mi,trial] = time() - tic;
    tic = time();
    Y   = S*X;
    Times[ALGO,Mi,trial] = time() - tic;

    ALGO = 2; # Fast JL, DCT
    tic     = time();
    D = spdiagm(0 => broadcast(sign,randn(M)) );
    ind     = rand(1:M, Int64(m));
    Times_setup[ALGO,Mi,trial] = time() - tic;
    tic = time();
    Y       = dct( D*X );
    Y       = Y[ind,:];
    Times[ALGO,Mi,trial] = time() - tic;

    ALGO = 3;  # Fast JL, Hadamard
    tic     = time();
    D = spdiagm(0 => broadcast(sign,randn(M)) );
    ind     = rand(1:M, Int64(m));
    Times_setup[ALGO,Mi,trial] = time() - tic;
    tic = time();
    Y       = M*fwht_natural( D*X, 1 );
    Y       = Y[ind,:];
    Times[ALGO,Mi,trial] = time() - tic;

    # ALGO = 4; # Count
    # tic     = time();
    # D = spdiagm(0 => broadcast(sign,randn(M)) );
    # indx_map        = Int64.(rand(1:m,M));
    # Times_setup[ALGO,Mi,trial] = time() - tic;
    # tic = time();
    # Y = slowCountSketch(X'*D,indx_map);
    # Times[ALGO,Mi,trial] = time() - tic;

    let ALGO = ALGO
    # Sparse. We can normalize later. Does that help speed?
      for Sparsity = 1:3
          ALGO        = ALGO + 1;
          s           = 2.0^(1-Sparsity)*sqrt(M);
          density     = 1/(2*s);
          ALGO_NAMES[ALGO] = string(round(100*density,digits=2), "% sparse");
          tic = time();
          S   = sprandn(Int64(m),M,density); # this takes longer than the multiply!
          S   = broadcast(sign, S);
          Times_setup[ALGO,Mi,trial]     = time() - tic;
          tic = time();
          Y   = sqrt(s)*(S*X);
          Times[ALGO,Mi,trial] = time() - tic;
      end # for loop for sparsity
    end # let block
  end # for loop for trials
end # for loop M list
# ---------------------------------------------------------------------------- #
# Plots
Data = Times;
mn  = reshape(mean(Data, dims = 3), nAlgos, length(M_list));
plot(M_list, mn', yscale = :log10, xscale = :log10, legend = :topleft,
        label = ALGO_NAMES, title = "times to apply sketch",
        titlefontsize = 10)
xlabel!("size M")
ylabel!("times in seconds")
y   = M_list/M_list[1];
plot!(M_list, mn[1].*y.^2, label = "M^2", linecolor = :black)
plot!(M_list, mn[3,1].*(M_list.*log.(M_list)/(M_list[1]*log.(M_list[1]))),
      label = "M log M", linecolor = :black, linestyle = :dash)
plot!(M_list, minimum(mn[:,1])*y, label = "M", linecolor = :black,
      linestyle = :dot)

# ---------------------------------------------------------------------------- #
Data = Times_setup + Times;
mn   = reshape(mean(Data,dims=3),nAlgos, length(M_list));;
plot(M_list, mn', yscale = :log10, xscale = :log10, legend = :topleft,
        label = ALGO_NAMES, title = "times to apply sketch with setup",
        titlefontsize = 10)
xlabel!("size M")
ylabel!("times in seconds")
y   = M_list/M_list[1];
plot!(M_list, mn[1].*y.^2, label = "M^2", linecolor = :black)
plot!(M_list, mn[3,1].*(M_list.*log.(M_list)/(M_list[1]*log.(M_list[1]))),
      label = "M log M", linecolor = :black, linestyle = :dash)
plot!(M_list, minimum(mn[:,1])*y, label = "M", linecolor = :black,
      linestyle = :dot)

# ---------------------------------------------------------------------------- #
