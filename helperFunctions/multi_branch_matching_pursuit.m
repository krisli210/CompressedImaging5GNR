function [x_est, err] = multi_branch_matching_pursuit(y, A, K, L, TOL)
% Inputs:
% y - measurement vector
% A - sensing matrix
% K - sparsity level of the input vector
% L - number of branches
% TOL - stopping criteria for residual error
%
% Outputs:
% x_est - estimated input vector
% err - final residual error

N = size(A, 2); % length of input vector

% Initialize variables
r = y; % residual
S = []; % support set
x_est = zeros(N, 1); % estimated input vector
err = norm(r); % initial residual error

for l = 1:L
    % Compute inner products
    inner_prods = abs(A'*r);
    
    % Sort the inner products in descending order
    [~, idx] = sort(inner_prods, 'descend');
    
    % Select the top K indices from the sorted list
    idx = idx(1:K);
    
    % Add selected indices to support set
    S = union(S, idx);
    
    % Solve least-squares problem to update estimated input vector
    x_S = A(:, S)\y;
    x_est(S) = x_S;
    
    % Update residual
    r = y - A*x_est;
    
    % Check stopping criteria
    err = norm(r);
    if err < TOL
        break;
    end
end

end
