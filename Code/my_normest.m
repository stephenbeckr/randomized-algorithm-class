function [e,cnt] = my_normest(S,St,n,tol, maxiter, nVectors)
%MY_NORMEST Estimate the matrix 2-norm. Improvement on matlab's version.
%   NORMEST(S) is an estimate of the 2-norm of the matrix S.
%   NORMEST(S,tol) uses relative error tol instead of 1.e-6.
%   [nrm,cnt] = NORMEST(..) also gives the number of iterations used.
%
%   NORMEST(S,St,n,[tol]) uses the funcions S() and St() to compute
%   the forward and transpose multiplication (where S is m x n)
%   This modification due to Stephen Becker, 10/01/09
%   If St is the empty matrix, then we assume S = S'
%   Note: this method is just the power method.
%
%   NORMEST( ..., maxiter ) uses at most maxiter iterations.
%   NORMEST( ..., maxiter, nVectors ) uses nVectors
%       and does the QR iteration
%
%   This function is intended primarily for sparse matrices,
%   although it works correctly and may be useful for large, full
%   matrices as well.  Use NORMEST when your problem is large
%   enough that NORM takes too long to compute and an approximate
%   norm is acceptable.
%
% Version modified by Stephen Becker to allow for function handles
%   Dec 2012, also suppors n=[n1,n2] for matrix-domain functions
% And modified Feb 2019 for vector-domain functions to do QR
%   iteration.

if nargin < 2, tol = 1.e-6; end
if nargin < 5, maxiter = 20; end
if nargin < 6, nVectors = []; end
IMPLICIT = false;
if isa(S,'function_handle')
    if isempty(St)
        St = S;  % we assume the matrix is symmetric;
    elseif ~isa(St,'function_handle')
        error('normest: must provide transpose function');
    end
    if nargin < 3
        error('normest: must provide width of matrix');
    end
    
    if nargin < 4, tol = 1.e-6; end
    IMPLICIT = true;
else
    if nargin >= 2 && isnumeric(St) && numel(St)==1, tol = St; end
    if nargin >= 3 && isnumeric(n) && numel(n)==1, maxiter = n; end
    n = size(S,2);
end

if isempty(nVectors)
    if numel(n) > 1
        nVectors = 1;
    else
        nVectors = 1;
    end
end

if ~IMPLICIT
    x = sum(abs(S),1)';
    if nVectors > 1
        x = [x, randn(n,nVectors-1)];
    end
else
%     x = ones(n,1); % can interact with some special operators
    if numel(n) == 1
        x = randn(n,nVectors);
    else
        if nVectors > 1
            error('Not compatible in this mode');
        end
        x = randn(n); % assume n is a size vector
    end
end

cnt = 0;
if nVectors > 1
    e = sqrt(max( sum(x.^2,1) ) );
    [x,~] = qr(x,0);
else
    e = norm(x(:));
    if e == 0, return, end
    x = x/e;
end
e0 = 0;
while abs(e-e0) > tol*e && cnt < maxiter
   e0 = e;
   if ~IMPLICIT
       Sx = S*x;
   else
       Sx = S(x);
   end
   if nnz(Sx) == 0
      Sx = rand(size(Sx));
   end
   if nVectors > 1
       e = sqrt(max( sum(Sx.^2,1) ) );
%        [Q,~] = qr(Sx,0);
%        Sx = Q;
   else
       e = norm(Sx(:));
   end
   if ~IMPLICIT
       x = S'*Sx;
   else
       x = St(Sx);
   end
   if nVectors > 1
      [Q,~] = qr(x,0);
      x     = Q;
   else
       x = x/norm(x(:));
   end
   cnt = cnt+1;
end
