function [MeanPattern, F_total] = genMeanPattern(U, NRB, N_T, azBins, BsArray, txCodebook)
    patterns = zeros(length(azBins), NRB * N_T);
    F_total = zeros(size(txCodebook, 1), 1);
    for ind = 1:NRB * N_T
        txAngles = randperm(size(txCodebook, 2), U);
        F = 1./sqrt(U) .* txCodebook(:, txAngles); % M x U
        F_total = F_total + sum(F, 2);
        patterns(:, ind) = pattern(BsArray, 28e9, azBins, Weights=sum(F, 2), Normalize=false, Type='power', EL=0);
    end
    MeanPattern = mean(patterns, 2);
end