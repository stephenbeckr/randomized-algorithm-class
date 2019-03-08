%{
Demonstration of the basic Randomized Kaczmarz Algorithm 
 (cf. "?A Randomized Kaczmarz Algorithm with Exponential Convergence"
   by Strohmer and Vershynin, 2008
   ?http://www.springerlink.com/index/10.1007/s00041-008-9030-4 )

For fancier versions, see:

- "?Acceleration of randomized Kaczmarz method via the
Johnson?Lindenstrauss Lemma" by Y. Eldar and D. Needell 2011
- "?Paved with good intentions: analysis of a randomized block kaczmarz
method" by D. Needell and J. Tropp 2012
- "?Stochastic gradient descent, weighted sampling, and the randomized
Kaczmarz algorithm" by D. Needell, N. Srebro and R. Ward 2016

Code: Stephen Becker 2019
%}


rng(0);
M   = 3e5;
N   = 1e2;
A   = randn(M,N);
x0  = randn(N,1);
b   = A*x0; % no noise, since we're not doing least squares, we're solving a system

tic
xLS     = A\b;
tm_LS = toc;
fprintf('Solving %d x %d system via Matlab classical method takes %.2f sec\n', ...
    M,N,tm_LS );

tic
[Q,R]   = qr(A,0);
xLS2    = R\(Q'*b);
tm_LS2  = toc;
fprintf('... or via a thin QR w/o pivoting takes %.2f sec\n', tm_LS2 );
%%
tic
rowNorms    = sum(A.^2,2);
% At = A'; % slow, but if I do this, then can accelerate iterations
tm_preprocess = toc;
% stem( rowNorms )
prob        = rowNorms/sum(rowNorms);
%%
x   = zeros(N,1);
maxIter     = 1e2;
errFcn      = @(x) norm(x - xLS );
errList     = zeros(maxIter,1);
tic
for k = 1:maxIter
%     i   = randsample(M,1,true,prob);
%     x   = x + (b(i)-A(i,:)*x)/rowNorms(i) * A(i,:)';
    
    iList = randsample(M,500,true,prob);
    for ind = 1:500
        i   = iList(ind);
        x   = x + (b(i)-A(i,:)*x)/rowNorms(i) * A(i,:)';
%         x   = x + (b(i)-At(:,i)'*x)/rowNorms(i) * At(:,i); % faster
    end
    
    errList(k)  = errFcn(x);
    if errList(k) < 1e-13
        break
    end
end
tm_Kaczmarz = toc;
%%
fprintf('Randomized Kaczmarz took %.2f sec = %.2f + %.2f sec; final error %.2e\n',...
    tm_preprocess + tm_Kaczmarz, tm_preprocess, tm_Kaczmarz, errList(k) );
%%
semilogy( errList,'o-','linewidth',2 )
xlabel('Epochs'); ylabel('Error'); set(gca,'fontsize',18); grid on
%% For a fair comparison, try with LSQR
tic
[xHat,flag,relres,iter]    = lsqr( A, b, 1e-13, 1e3 );
tm_CG = toc;
fprintf('LSQR took %.2f sec in %d iterations; final error %.2e\n',...
    tm_CG, iter, errFcn(xHat) );