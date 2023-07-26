function W = optEladLS(D, P_init, X, N)

    %% Elad Algorithm
    % Fixed parameters taken from Elad Fig. 3 - open to optimization
    t = .2; % Thresholding 
    gamma = .8; % Scaling
    p = size(P_init, 1); % number of measurements
    
    Q = 100; % Number of iterations
    coherences = zeros(1, Q);
    coherence_init = mean_coherence(P_init*D)*ones(1,Q);
    
    lambda = 10; % ridge regression scaling
    P_q = P_init;
    % P_q = randn(size(P_init));
    for q = 1:Q
        D_hat = P_q*D;
        D_hat = D_hat./sqrt(sum(abs(D_hat).^2)); % Norm columns to 1 - correct normalization??
        
        Gram = D_hat' * D_hat;
        t = abs(findThreshold(Gram, .10));

        %Shrinking
        mask1 = abs(Gram) >= t & ~eye(size(Gram, 1));
        % mask2 = abs(Gram) < t & abs(Gram) >= gamma*t;
        % mask3 = abs(Gram) < gamma*t;

        Gram_hat = Gram;
        Gram_hat(mask1) = gamma*Gram(mask1);
        % Gram_hat(mask2) = gamma*t * sign(Gram(mask2));
        % Gram_hat(mask3) = Gram(mask3);
        % Gram_hat = Gram_hat - diag(diag(Gram_hat)) + diag(diag(Gram)); % Replace the diagonal elements 

        % Step 5+6 
        S_q = lowRankSqrt(Gram_hat, p);
        
        % Step 6 find S_q s.t. S_q.'*S_q = Gram_hat
        
        % Step 7 Find P_{q+1} as min_P ||S_q - PD||^2_F
        P_q = S_q*pinv(D); % this isn't necessarilly supposed to yield low P_q-Pinit error

        coherences(q) = mean_coherence(P_q*D);
    end

    %% Determined P, now must solve LS to acquire W from P = kron(X^T, W)
    % Given in F. Roemer's PhD thesis sec 3.4
    % I know this implemented correctly because we can set Q = 0 and yield
    % W back

    M = size(X, 1);
    Ns = size(X, 2);
    X_block = zeros(M*N*Ns, N);
    I_N = eye(N);
    P_q = P_q./sqrt(sum(abs(P_q).^2));
    for m = 1:M
        start_ind = (m-1)*N*Ns+1;
        end_ind = start_ind + N*Ns - 1;
        X_block(start_ind:end_ind, :) = kron(I_N, X(m, :).');
    end
    
    % vec_W = kron(pinv(X_block), I_N) * reshape(P_q, [numel(P_q), 1]);

    % vec_W = norm(X, "fro")^-2 * kron(X_block', I_N) * reshape(P_q, [numel(P_q), 1]); %equivalent
    
    % X_block_block = kron(X_block, I_N);
    vec_W = kron((X_block'*X_block + lambda*I_N)^-1 * X_block', I_N) * reshape(P_q, [numel(P_q), 1]);

    W = reshape(vec_W, [N N]);

    % W = W + diag(max(W, [], 'all').*ones(1, N));
    % W = W./sqrt(sum(abs(W).^2));
end

function t = findThreshold(Gram, percentile)
    offDiag = Gram - diag(diag(Gram));
    A = sort(reshape(offDiag, [numel(offDiag), 1]),'descend');
    
    ind = floor(length(A)*percentile);
    t = A(ind);
end

function S = lowRankSqrt(G, p)
    [U, S, V] = svd(G);

    Up = U(1:p, 1:p);
    Sp = sqrt(S(1:p, 1:p));
    Vp = V(:, 1:p);

    S = Up * Sp * Vp';
end
