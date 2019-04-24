function C  = AMS_sketch( X, w, d, varargin )
% C  = AMS_sketch( X, w, d )
%   returns a linear sketch of the input matrix X, X is p x n
%       (note convention: *columns* not *rows* of X are the data pts)
%       and the output sketch C is of size w*d x n
%       where C(:,i)    = AMS_sketch( X(:,i) )
%
%   w controls number of buckets per hash (e.g., 2^8, 2^12)
%   d controls how many hashes we do (e.g., 7, or ceil(log2(1/.01)) )
%
% C  = AMS_sketch( X, w, d, parameters, values )
%   gives more options, e.g.,
%   
%   'saltSeed', seed
%   gives a new seed to the random number generator (default: 0)
%   which controls both hashes
%
%   'transposedX', false
%   if True, assumes X is n x p not p x n
%
% Stephen Becker
%   This version is good for matrices that are not too sparse,
%   and have more than a few columns n. If you are applying
%   this to a very sparse vector, then you really should
%   use hash functions if you want sub-linear time
%   (this implementation calls several Count sketches, 
%    which I have implemented in a way that is not sub-linear
%    for a single column).

prs = inputParser;
addParameter(prs,'saltSeed',.0);
addParameter(prs,'transposedX',false);
parse(prs,varargin{:});
saltSeed        = prs.Results.saltSeed;
transposedX     = prs.Results.transposedX;


if transposedX
    [n,p]   = size(X);
else
    [p,n]   = size(X);
end
% if p > intmax('uint32')
%     error('Dimensions of input matrix are too large!');
% end
% if w > intmax('uint16')
%     error('Code needs to be updated if you want w > 2^16');
% end
% if 2*d > 20
%     error('Code needs to be updated if you want 2*d > 20');
% end

rng( saltSeed );
% saltPerm    = randperm( 20 ); % output of SHA has 20 int8's
% 
% 
% 
% % C   = zeros( d, w, n );
% C   = zeros( d*w, n );
% 
% for j = 1:p 
%     Engine  = java.security.MessageDigest.getInstance('SHA');
%     Engine.update(typecast( uint32(j), 'uint8'));
%     L       = uint16(typecast( Engine.digest, 'uint8' ));
%     L       = L( saltPerm );
%     binaryHashes    = sign( randn(d,1) ); % don't even need a hash "function"
%     for k = 1:d
%         ell     = L(2*k-1)*2^8 + L(2*k);
%         ell     = mod( ell, w ) + 1;    % keep it in range, and make it 0-based
%         ind     = sub2ind( [d,w], k, ell );
%         C ( ind, : ) = C( ind, : ) + binaryHashes(k)*X( j, : );
%         %C ( k, ell, : ) = C( k, ell, : ) + X( j, : ); % same idea, if C is
%         %   a tensor
%     end
% end
% 

% % Make it faster... by skipping the hash!
% Instead, just call Count Sketch function


m       = w;
M       = p;
useTranspose    = true;
C       = zeros( d*w, n );
for k = 1:d
    D        = spdiags(sign(randn(M,1)),0,M,M); % bsxfun() is another efficient way to do this
    indx_map = int64(randi(m,M,1));
    if transposedX
        % I want X', but it's already transposed, so don't worry!
        C( (1+(k-1)*w):k*w, :) =  countSketch_BLAS(X*D,indx_map,m,useTranspose)';
    else
        C( (1+(k-1)*w):k*w, :) =  countSketch_BLAS(X'*D,indx_map,m,useTranspose)';
    end
end