function out = kr(A,B)
% Implements the Khatri-Rao product between two matrices A and B
% Input: A - m x n matrix
%        B - p x n matrix
% Output: out - mp x n matrix

% Get the dimensions of the matrices
[m, n] = size(A);
[p, ~] = size(B);

% Initialize the output matrix
out = zeros(m*p, n);

% Compute the Khatri-Rao product
for j = 1:n
    out(:,j) = kron(A(:,j), B(:,j));
end
end
