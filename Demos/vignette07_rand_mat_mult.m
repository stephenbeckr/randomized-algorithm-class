% function vignette07_rand_mat_mult
%{
Vignette #7
APPM 4720/5720 Randomized Algorithms, Spring 2019

Demo to show how to approximate a matrix-matrix product using randomization
based on the approach of P. Drineas and M. Mahoney. See section four of
"Lectures on Randomized Linear Algebra" by Drineas and Mahoney for additional
details and analysis (?http://arxiv.org/abs/1712.08880)

Algorithm randomized matrix multiplication.
Inputs:
 A is m x n
 B is n x p
 integer  1 <= c <= n
 probability distribution of length n

repeat for k = 1,\dots, c:
  1. pick i \in {1,\dots,n} with P(i = k) = p_k iid with replacement.
  2. set C(:,k) = 1/sqrt(c * p_i)*A(:,i) and R(k,:) = 1/sqrt(c * p_i)*B(i,:).

return C, R, and CR = sum_{k=1}^c 1/(c * p_i) * A(:,i) B(i,:).

Code by Jake Knigge, modifications by Stephen Becker

Exercise: can you use leverage scores to improve sampling?
  Issue #1: you have to balance leverage scores of A *and* of B
  Issue #2: what happens if n > p or n > m?
%}
% ---------------------------------------------------------------------------- %
rng(2);                                     % set seed for reproducibility
n_sims = 1000; m = 500; n = 20; p = 250;   % parameters for matrix sizes
c = min(max(ceil(0.5*n),1),n);              % subsampling parameter
% storage for simulations
fro_norms = zeros(n_sims, 1); fro_norms_opt = fro_norms;
% ---------------------------------------------------------------------------- %
% generate "data" matrices
A = randn(m,n)/sqrt(n); A(:,1) = 5 + 5*rand(m,1); % make the first column "big"
B = randn(n,p)/sqrt(n); AB = A*B;
% ell_2 norms: columns of A and rows of B
col_norm_A = sqrt(sum(A.^2, 1)); row_norm_B = sqrt(sum(B.^2, 2));
% MATLAB R2018b includes a 'vecnorm' function to calculate column- or row-wise
% norms of matrices; see https://www.mathworks.com/help/matlab/ref/vecnorm.html.
% ---------------------------------------------------------------------------- %
% probabilities for sampling---naive = uniform; optimal uses norm information
probs = ones(n, 1)/n;
probs_opt = ( (col_norm_A .* row_norm_B') / (col_norm_A * row_norm_B) )';
% compute theoretical upper bounds on expected squared Frobenius norms
upper_bound = sum(col_norm_A.^2 .* (row_norm_B.^2)' ./ (c*probs'));
upper_bound_opt = (col_norm_A * row_norm_B)^2 / c;
% ---------------------------------------------------------------------------- %
% simulation
for t = 1:n_sims
    % take samples
    samplesUniform  = randperm( n, c );
    replace         = true;
    samplesOptimal  = randsample( n, c, replace, probs_opt );
    rescale         = 1./sqrt( c/n );
    rescale_opt     = 1./sqrt( c*probs_opt( samplesOptimal) );
    C               = A(:,samplesUniform)*diag( rescale );
    C_opt           = A(:,samplesOptimal)*diag( rescale_opt );
    R               = diag( rescale )    *B(samplesUniform,:);
    R_opt           = diag( rescale_opt )*B(samplesOptimal,:);
    
    % compute random matrix product via outerproduct
    CR      = C*R; 
    CR_opt  = C_opt*R_opt;
    % calculate Frobenius norms of actual vs. randomized
    fro_norms(t)    = norm(AB-CR,'fro'); 
    fro_norms_opt(t)= norm(AB-CR_opt,'fro');
end 

% display comparisons of simulations
formatSpec = '||AB-CR||_F simulation %4.2f vs. upper bound %4.2f (naive)\n';
fprintf(formatSpec, sqrt(mean(fro_norms.^2)), sqrt(upper_bound))
formatSpec = '||AB-CR||_F simulation %4.2f vs. upper bound %4.2f (optimal)\n';
fprintf(formatSpec, sqrt(mean(fro_norms_opt.^2)), sqrt(upper_bound_opt))

% fprintf('relative err: %.3f (naive), %.3f (optimal)\n', ...
    

% NOTE: the Python verison of this script compares the squared Frobenius norms
% from the simulations vs. the upper bounds. The comparisons in this script are
% different (and make use Jensen's inequality for concave functions), but their
% interpretations are the same.
% end % vignette_rand_mat_mult function
% ---------------------------------------------------------------------------- %
% If you don't have randsample (in new versions of Stat toolbox)
%   then repeatedly call this randSample function
% helper function
% function y = randSample(x)
% % function randSample
% % samples an integer from a prob. distribution x (i.e., x >= 0, sum(x) = 1)
% idx = (1:length(x))'; cdf = cumsum(x); z = rand(); y = min(idx(z <= cdf));
% end
% ---------------------------------------------------------------------------- %
% test of randSample() function using optimal sampling probabilities from the
% above code
% y = zeros(1e5,1);
% for t = 1:(length(y))
%     y(t) = randSample(probs_opt);
% end
% histogram(y,'Normalization','probability'); hold on, plot(1:n,probs_opt,'*r');
