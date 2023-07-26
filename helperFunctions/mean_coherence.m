function mu = mean_coherence(A)
    
    [m, n] = size(A);

    nn = sqrt(sum(A.*conj(A),1));
    nA = bsxfun(@rdivide,A,nn);  % nA is a matrix with normalized columns

    innerProducts = abs(nA' * nA);
    innerProducts(logical(eye(n))) = 0;  % Set diagonal elements to zero
    
    mu = sum(innerProducts(:)) / (n * (n - 1));
end