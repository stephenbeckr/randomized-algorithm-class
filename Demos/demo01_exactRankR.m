%{

Vignette #1
APPM 4720/5720 Randomized Algorithms, Spring 2019
Stephen Becker

This is not a practical algorithm, since it only works
if the matrix is *exactly* low rank

%}

rng(0);  % seed the random number generator so this is reproducible

% -- Generate a low-rank matrix (m x n, with rank r )
n   = 4e3;
m   = n;
r   = 100;  % rank

A   = randn(m,r)*randn(r,n);

%% -- Find its SVD with conventional methods
tic
[U,Sigma,V]     = svd(A,'econ');
toc
% 35 seconds
% semilogy( diag(Sigma), 'o-' )
% Matlab doesn't understand that the matrix is not full rank

clear U Sigma V
%% -- Find its SVD with a randomized method
tt  = tic;
tic; Omega   = randn(n,r);   toc
tic; Y       = A*Omega;      toc
tic; [Q,R]   = qr(Y,0);      toc
tic; QtA     = Q'*A;         toc
tm  = toc(tt);

A_estimate  = Q*QtA;
err         = norm( A - A_estimate, 'fro' )/norm(A,'fro');
fprintf('||A-A_estimate||_F/||A||_F is %g\n', err );
fprintf('Overall time: %g seconds\n', tm );