function W_q = optEladLS_v2(D, P_init, X, N)

    %% Elad Algorithm
    % Fixed parameters taken from Elad Fig. 3 - open to optimization
    t = .1; % Thresholding 
    gamma = .8; % Scaling
    p = size(P_init, 1); % number of measurements
    
    Q = 50; % Number of iterations
    coherences = zeros(1, Q);

    % W_q = eye(N);
    W_q = randn(N)+1j*randn(N);
    W_q = W_q./sqrt(sum(abs(W_q).^2));

    M = size(X, 1);
    Ns = size(X, 2);
    X_block = zeros(M*N*Ns, N);
    I_N = eye(N);
    for m = 1:M
        start_ind = (m-1)*N*Ns+1;
        end_ind = start_ind + N*Ns - 1;
        X_block(start_ind:end_ind, :) = kron(I_N, X(m, :).');
    end

    % P_q = randn(size(P_init));
    for q = 1:Q

        P_q = kron(X.', W_q);
        coherences(q) = mean_coherence(P_q*D);

        D_hat = P_q*D;
        % D_hat = D_hat./sqrt(sum(abs(D_hat).^2)); % Norm columns to 1 - correct normalization??
        
        Gram = D_hat' * D_hat;
        t = abs(findThreshold(Gram, .10));

        %Shrinking
        mask1 = abs(Gram) >= t & ~eye(size(Gram, 1));
        % mask1 = mask1 & ( rand(size(mask1)) > .5);

        Gram_hat = Gram;
        Gram_hat(mask1) = gamma*Gram(mask1);


        % Step 5+6 
        S_q = lowRankSqrt(Gram_hat, p);
        % norm(S_q'*S_q - Gram_hat, 'fro')
       
        % Step 6 find S_q s.t. S_q'*S_q = Gram_hat
        
        % Step 7 Find P_{q+1} as min_P ||S_q - PD||^2_F
        P_q = S_q*pinv(D); % % The error here lies in cutting off the last p > k entries. You can use imagesc() to see this. 
        % norm(S_q - P_q*D, 'fro')
        % norm(kron(X.',W_q)-P_q, 'fro')

        vec_W_q = norm(X, "fro")^-2 * kron(X_block', I_N) * reshape(P_q, [numel(P_q), 1]); %equivalent
        W_q = reshape(vec_W_q, [N N]).';
        W_q = W_q./sqrt(sum(abs(W_q).^2));
    end

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
