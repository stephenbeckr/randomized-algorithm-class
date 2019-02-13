function y = Hadamard_teaching_code(x)
% y = Hadamard_teaching_code(x)
%   applies the Hadamard transform to x
%   If x has more than one column, the transform is applied
%   to each column.
% This code is not fast, but it shows you how to exploit
%   the structure of the transform.
% Note: this code does not do any sub-sampling

[m,n]   = size(x);
if 2^nextpow2(m) ~= m
    error('Must have leading dimension of x be power of 2');
end

y   = x;
for bit = 1:log2(m)
    k   = 2^bit; % e.g., 2, 4, ..., m
    k2  = 2^(bit-1); % e.g., 1, 2, ..., m/2
    
    y   = reshape( y, k, [], n );
    tmp = y(1:k2,:,:);
    y(1:k2,:,:)     = y(1:k2,:,:) + y(k2+1:k,:,:);
    y(k2+1:k,:,:)   = tmp         - y(k2+1:k,:,:);
    y   = reshape( y, m, n);
end