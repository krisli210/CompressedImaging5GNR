function x_hat = omp(A, y, max_iter, residual_thresh)

% A: Sensing matrix of size M x N
% y: Measurement vector of size M x 1
% K: Sparsity level
% max_iter: Maximum number of iterations
% residual_thresh: Threshold residual

[M, N] = size(A);
S = []; % Support set
r = y; % Residual
x_hat = zeros(N, 1); % Estimated signal
iter = 0; % Iteration counter

while norm(r) > residual_thresh && iter < max_iter
    % Select atom that maximizes correlation with residual
    [~, idx] = max(abs(A' * r));
    
    % Add selected index to support set
    S = [S, idx];
    
    % Solve least-squares problem over support set
    A_S = A(:, S);
    x_S = A_S \ y;
    
    % Update estimated signal and residual
    x_hat(S) = x_S;
    r = y - A * x_hat;
    
    % Increment iteration counter
    iter = iter + 1;
end
