%{
Demo of the Iterative Hessian Sketch (IHS) cf. Pilanci and Wainwright
and of the preconditioned approaches (BLENDENPIK, LSRN)

These are two methods to get high-accuracy l2 regression

The goal is to approximate the solution of
  min_{x} || Ax-b ||^2
where
  A is M x N
and we are assuming M >> N.

%}

addpath ~/Repos/randomized-algorithm-class/Code/
rng(0);

M   = 8e4;
N   = 5e2;
% M   = 5e3; N = 20;
A   = randn(M,N)*diag(logspace(0,3,N))*(randn(N)+.1*eye(N));
fprintf('Condition number is %.3e\n', cond(A)  )

x   = randn(N,1);
b   = A*x;
b   = b + .3*norm(b)/sqrt(M)*randn(M,1);

%% Solve via dense method, about 7 seconds
tic
xLS     = A\b;
tm_LS = toc;
fprintf('\nSolved via classical least-squares in %.2f seconds\n\n',tm_LS);

%% Solve via another dense method, not as safe if ill=conditioned
tic
[Q,R]   = qr(A,0);
xHat    = R\(Q'*b);
tm_LS_2 = toc;

err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
err2        = norm( xHat - xLS)/norm(xLS);
err3        = norm( A*(xHat - xLS))/norm(A*xLS);

fprintf('== Classical algorithm 1, QR without pivoting, less robust but faster ==\n');
fprintf('Took %.2f sec, err metrics %.2e and %.2e and %.2e\n\n', ...
    tm_LS_2, err1,err2,err3);

% At least do pivoting
tic
[Q,R,e]     = qr(A,0);
xHat(e)       = R\(Q'*b); %
tm_LS_3 = toc;

err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
err2        = norm( xHat - xLS)/norm(xLS);
err3        = norm( A*(xHat - xLS))/norm(A*xLS);

fprintf('== Classical algorithm 2, QR with column pivoting ==\n');
fprintf('Took %.2f sec, err metrics %.2e and %.2e and %.2e\n\n', ...
    tm_LS_3, err1,err2,err3);

%% Take a sketch
rng(2);
m   = 40*N;
sketchType = {'count','FJLT'}; % Gaussian is too slow!
for i = 1:2
    type    = sketchType{i};
    tic;
    sketchFcn   = sketch( m, M, type );
    SAb         = sketchFcn([A,b]);
    time_preprocess     = toc;
    fprintf(' -- Sketch type %s took %.2f sec\n', type, time_preprocess );
end


%% Try a normal style sketches

tic;
SA          = SAb(:,1:N);
Sb          = SAb(:,N+1);
xHat        = SA\Sb;
time_solve          = toc;
err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
err2        = norm( xHat - xLS)/norm(xLS);
err3        = norm( A*(xHat - xLS))/norm(A*xLS);
fprintf('== Standard JL style sketch ==\n');
fprintf('Took %.2f = %.2f + %.2f sec, err metrics %.2e and %.2e and %.2e\n\n', ...
    time_preprocess+time_solve, time_preprocess,time_solve, err1,err2,err3);

%% Try the iterative Hessian Sketch (run the above block first to get SAb)
nBlocks     = 4;
mm          = floor(m/nBlocks);
tic
xHat        = zeros(N,1);
bHat        = b;
fprintf('== Iterative Hessian Sketch ==\n');
tic
for i = 1:nBlocks
    startInd    = 1 + (i-1)*mm;
    endInd      = i*mm;
    SA          = sqrt(m/mm)*SAb(startInd:endInd,1:N); % renormalize!
%     Sb          = SAb(startInd:endInd,N+1);
%     xx          = SA\Sb;  $ regular sketching
    xx          = (SA'*SA )\(A'*bHat);
    xHat        = xHat + xx;
    bHat        = bHat - A*xx;
    err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
    err2        = norm( xHat - xLS)/norm(xLS);
    err3        = norm( A*(xHat - xLS))/norm(A*xLS); % need < 1
    fprintf('  contraction factor at iter %d is %.4f\n', i, err3 );
end
time_solve_IHS = toc;
fprintf('Took %.2f = %.2f + %.2f sec, err metrics %.2e and %.2e and %.2e\n\n', ...
    time_preprocess+time_solve_IHS, time_preprocess,time_solve_IHS, err1,err2,err3);

%% Try preconditioning
fprintf('== Computing thin QR on sketched data ==\n');
tic
[Q,R]       = qr( SA, 0 ); % thin QR
time_QR     = toc;

k1=cond( SA/R ); % unless we had precision issues, this ought to be 1
k2=cond( A/R );  % and this thing we *hope* is small
% Note: cond( A/R ) is a nicer way to write cond( A*inv(R) )

fprintf('QR on SA took %.2f sec, cond(SA*inv(R)) is %.2f, cond(A*inv(R)) is %.2f\n\n',...
    time_QR,k1,k2);
%% For reference, use LSQR to solve, without preconditioning
fprintf('== LSRN, for reference, with 100 iterations ==\n');
tol     = 1e-8;
maxit   = 1e2;
tic
xHat      = lsqr(A,b,tol,maxit);
time_LSQR = toc;
err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
err2        = norm( xHat - xLS)/norm(xLS);
err3        = norm( A*(xHat - xLS))/norm(A*xLS);
fprintf('Took %.2f sec for LSQR, err metrics %.2e and %.2e and %.2e\n\n', ...
    time_LSQR, err1,err2,err3);


fprintf('== LSRN, for reference, with 500 iterations ==\n');
tol     = 1e-8;
maxit   = 5e2;
tic
xHat      = lsqr(A,b,tol,maxit);
time_LSQR = toc;
err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
err2        = norm( xHat - xLS)/norm(xLS);
err3        = norm( A*(xHat - xLS))/norm(A*xLS);
fprintf('Took %.2f sec for LSQR, err metrics %.2e and %.2e and %.2e\n\n', ...
    time_LSQR, err1,err2,err3);

%% Now try LSQR with preconditioning
fprintf('== Preconditioned LSQR (a la BLENDENPIK/LSRN) ==\n');
tol     = 1e-9;
maxit   = 1e2;
tic
xHat      = lsqr(A,b,tol,maxit,R);
time_LSQR_R = toc;
err1        = norm( A*xHat - b )/norm(A*xLS-b) - 1;
err2        = norm( xHat - xLS)/norm(xLS);
err3        = norm( A*(xHat - xLS))/norm(A*xLS);
fprintf('Took %.2f = %.2f + %.2f + %.2f sec for LSQR, err metrics %.2e and %.2e and %.2e\n\n', ...
    time_preprocess+time_QR+time_LSQR_R, time_preprocess,time_QR,time_LSQR_R, err1,err2,err3);
