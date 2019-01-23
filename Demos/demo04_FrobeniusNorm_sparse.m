%{

Accessing rows vs accessing columns of a sparse matrix

We saw in demo03 that for a dense matrix, in one particular
code and size of matrix, there was about a factor of 3 or 4x
improvement if we looped over columns vs rows

What about for sparse matrices?

%}

%% Generate a sparse matrix with random entries, 20% sparse
A = sprandn(1e4,1e4,0.2); % this is 10,000 x 10,000

%% 
disp('Compute Frobenius norm by looping over columns');
s=0; 
tic
for i=1:size(A,2)
    s = s + norm( A(:,i) )^2; 
end
t1=toc;
fprintf('Frobenius norm is %e, took %f seconds\n\n', s^2, t1 );

%% 
disp('Compute Frobenius norm by looping over rows')
s=0; 
tic
for i=1:size(A,1)
    s = s + norm( A(i,:) )^2; 
end
t2 = toc;
fprintf('Frobenius norm is %e, took %f seconds\n', s^2, t2 );

fprintf('Access via column vs row is %.1fx faster\n', t2/t1 );


%% Extra: can you tell what's the difference?
% Why do we get different values? Why is one slower?
X   = randn(5e3);
tic; s1=norm(X(:)); toc
tic; s2=norm(X); toc