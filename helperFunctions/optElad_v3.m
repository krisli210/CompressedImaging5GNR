function P_q = optElad_v3(D, P_init)

    %% Elad Algorithm
    % Fixed parameters taken from Elad Fig. 3 - open to optimization
    % Taken from now optimizing WH_RX
    t = .2; % Thresholding 
    gamma = .8; % Scaling
    p = size(P_init, 1); % number of measurements
    
    Q = 50; % Number of iterations
    coherences = zeros(1, Q);
    coherence_init = mean_coherence(P_init*D)*ones(1,Q);
    
    P_q = P_init;
    % P_q = randn(size(P_init));
    for q = 1:Q
        D_hat = P_q*D;
        % D_hat = D_hat./sqrt(sum(abs(D_hat).^2)); % Norm columns to 1 - correct normalization??
        
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
