function vignette_rsvd
%{
Vignette - randomized singular value decomposition (RSVD).
APPM 4720/5720 Randomized Algorithms, Spring 2019

Demo to illustrate a simple RSVD based on the two-stage approach of Martinsson. 
In particular, see section four of "Randomized Methods for Matrix Computations"
by P.G. Martinsson as appearing in "The Mathematics of Data" for details and 
additional analysis. An updated version of the survey paper can be found at:
https://arxiv.org/pdf/1607.01649.pdf. (Again, see section four, pages 8-9.)

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
%}
% ---------------------------------------------------------------------------- %
rng(2);                             % set seed for reproducibility
n_sims = 5;                         % number of simulations
n_subs = 25;                        % number of subsamples
m = 2000;                           % rows of matrix
% n_sims = 50; n_subs = 25; m = 1000;         % alt. run parameters
% n_sims = 1;  n_subs = 25; m = 5000;         % alt. run parameters
n = 15*ceil(log(m));                % columns of matrix
k = ceil(log(m));                   % target rank
p = max([ceil(log(m)) 10]);         % oversampling parameter
% ---------------------------------------------------------------------------- %
times = zeros(n_sims, 2);           % bookkeeping for run times
fro_mean = zeros(n_sims, n_subs);   % bookkeeping for mean Frobenius norm
op_mean = zeros(n_sims, n_subs);    % bookkeeping for mean operator norm
norm_bound = zeros(n_sims, 2);      % bookkeeping for theoretrical norm
fro_lo = zeros(n_sims, 1);          % bookkeeping for Frobenius norm bound
op_lo = zeros(n_sims, 1);           % bookkeeping for operator norm bound
% ---------------------------------------------------------------------------- %
for j = 1:n_sims
    A = [2+rand(m, 2) randn(m, k-2) 0.01*randn(m, n-k)]/sqrt(m); % "data" matrix
    % NOTE: the matrix A has k columns that lead to "important" singular values.
    % The remaing n-k columns correspond to fast-decaying singular values.
    for i = 1:n_subs
        % stage A of randomized SVD
        G = randn(n, k+p);                              % Gaussian random matrix
        Y = A*G;                                        % sample matrix
        [Q ~] = qr(Y,0);                                % orthonormalize Y
        % stage B of randomized SVD
        B = Q'*A;                                       % form small matrix
        tic; [U_B, D_B, V_B] = svd(B, 'econ'); times(j, 1) = toc;   % SVD of B
        U = Q*U_B;                                      % rank k matrix
        % bookkeeping and comparisons
        tic; [U_A, D_A, V_A] = svd(A, 'econ'); times(j, 2) = toc;   % SVD of A
        fro_mean(j,i) = norm(A - Q*Q'*A, 'fro');
        op_mean(j,i) = norm(A - Q*Q'*A);
    end % end of inner simulation loop (i.e., simulation for a fixed matrix A)
    fro_lo(j, 1) = sum(diag(D_A(k+1:min([m,n]),k+1:min([m,n]))).^2)^0.5;
    op_lo(j, 1) =  D_A(k+1,k+1);
    norm_bound(j, 1) = (1 + k / (p-1))^(0.5) * fro_lo(j, 1);     % Frobenius 
    norm_bound(j, 2) = (1 + sqrt(k / (p-1))) * op_lo(j, 1) + ... % operator
                            exp(1) * (sqrt(k+p) / p) * fro_lo(j, 1);
end % end of outer simulation loop
% ---------------------------------------------------------------------------- %
figure; % plot showing singular value decay
semilogy(diag(D_A), '-bo'), hold on, semilogy(diag(D_B), '--xr'), hold off;
title('singular value comparison/decay (final simulation)'); 
legend('full','randomized')
% ---------------------------------------------------------------------------- %
figure; subplot(2,1,1); % plots comparing computation times for SVDs
semilogy(times(:,1), '--xr'), hold on, semilogy(times(:,2), '-bo'), hold off;
legend('randomized', 'full', 'Location', 'east')
title('SVD timing comparison')
subplot(2,1,2); 
plot(times(:,2) ./ times(:,1), '-ks');
title('ratio of SVD times'); 
legend('full-to-randomized', 'Location', 'southeast');
% ---------------------------------------------------------------------------- %
figure; % plots of average Frobenius and operator norms vs. theoretical bounds
subplot(1,2,1);
plot(norm_bound(:,1),'-ob'), hold on, plot(mean(fro_mean, 2),'--xr');
        plot(fro_lo, ':*b'); ylim([0 max(max(norm_bound))*1.1]); hold off; 
        title('Eckhart-Young bounds: Frobenius')
        legend('upper','mean','lower','Location','southoutside')
subplot(1,2,2);
plot(norm_bound(:,2),'-ob'), hold on, plot(mean(op_mean, 2),'--xr');
        plot(op_lo, ':*b'); ylim([0 max(max(norm_bound))*1.1]); hold off; 
        title('E-Y bounds: operator')
        legend('upper','mean','lower','Location','southoutside')
end % end of function
% ---------------------------------------------------------------------------- %