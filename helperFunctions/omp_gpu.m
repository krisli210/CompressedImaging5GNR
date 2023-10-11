function x_hat = omp_gpu(A, y, max_iter, residual_thresh)

% A: Sensing matrix of size M x N
% y: Measurement vector of size M x 1
% max_iter: Maximum number of iterations
% residual_thresh: Threshold residual

[M, N] = size(A);
S = gpuArray([]); % Support set on GPU
r = gpuArray(y); % Residual on GPU
A = gpuArray(A);
x_hat = gpuArray.zeros(N, 1); % Estimated signal on GPU
iter = 0; % Iteration counter

while norm(r) > residual_thresh && iter < max_iter
    % Select atom that maximizes correlation with residual
    [~, idx] = max(abs(A' * r));
    
    % Add selected index to support set if it's not already in there
    if ~ismember(idx, S)
        S = [S, idx];
    end
    
    % Solve least-squares problem over support set
    A_S = A(:, S);
    x_S = A_S \ y;
    
    % Update estimated signal and residual
    x_hat(S) = x_S;
    r = y - A * x_hat;
    
    % Increment iteration counter
    iter = iter + 1;
end

% Transfer the result back to CPU
x_hat = gather(x_hat);
