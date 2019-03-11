function [D,I]   = pdist2_faster(X,Y,style,smallStr,K)
%  pdist2_faster Pairwise distance between two sets of observations.
%     D = pdist2(X,Y) returns a matrix D containing the Euclidean distances
%     between each pair of observations in the MX-by-N data matrix X and
%     MY-by-N data matrix Y. Rows of X and Y correspond to observations,
%     and columns correspond to variables. D is an MX-by-MY matrix, with the
%     (I,J) entry equal to distance between observation I in X and
%     observation J in Y.
%  
%     D = pdist2(X,Y,DISTANCE) computes D using DISTANCE.  Choices are:
%  
%         'euclidean'        - Euclidean distance (default)
%         'squaredeuclidean' - Squared Euclidean distance
%
% Also gives the index corresponding to the smallest entry if requested ...
%   see documentation for pdist2 from the Statistics Toolbox
% Code by Stephen Becker, March 2019
% See also pdist2

if nargin < 2 || isempty(Y)
    Y   = X;
end

[Mx,N]  = size(X);
[My,N]  = size(Y);

XtY     = X*Y';
nrm1    = sum(X.^2,2);
nrm2    = sum(Y.^2,2);

D       = nrm1*ones(1,My) + ones(Mx,1)*nrm2' - 2*XtY;

if nargin > 2 && ~isempty(style)
    if isa(style,'function_handle')
        % apply the function handle!
        D   = style(D);
    else
        switch lower(style)
            case 'euclidean'
                D   = sqrt(D);
            case 'squaredeuclidean'
        end
    end
else
    % by default, Euclidean
    D   = sqrt(D);
end

% and similar to pdist2,
if nargin > 3
    if nargin ==5 && ~isempty(K)
        if K~=1
            error('K ~= 1 not supported');
        end
    end
    if strcmpi(smallStr,'smallest')
        [~,I]     = min( D, [], 1 ); % over 1st dimension, like Matlab's convention
        D         = D(I);
    else
        fprintf('Bad string: %s\n');
        error('Wrong 4th input');
    end
end
end