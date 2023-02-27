function [mag_est, I] = solveCS_OMP(y, Phi, Psi)
    
    %returns estimated magnitudes mag_est under support I 
    A = Phi * Psi;
    sensingDict = sensingDictionary('CustomDictionary', A);
    [mag_est, YI, I, R] = matchingPursuit(sensingDict, y, maxIterations=10, Algorithm="OMP", maxerr={"L1", 1e-4});
    
    A_sub = Phi * Psi(:, I); % Solve magnitude posthence via direct linsolve against estimated support
    mags = linsolve(A_sub, y);
    mag_est(I) = mags;
end