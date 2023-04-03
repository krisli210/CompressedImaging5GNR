function [mag_est, I, YI, R] = solveCS_OMP(y, Phi, Psi, varargin)
    if (nargin == 4)
        OMP_iters = varargin{1};
    else
        OMP_iters = 10;
    end
    
    %returns estimated magnitudes mag_est under support I 
    A = Phi * Psi;
    sensingDict = sensingDictionary('CustomDictionary', A);
    [mag_est, YI, I, R] = matchingPursuit(sensingDict, y, maxIterations=OMP_iters, Algorithm="OMP", maxerr={"L1", 1e-15});
    
    A_sub = Phi * Psi(:, I); % Solve magnitude posthence via direct linsolve against estimated support
    mags = linsolve(A_sub, y);

    mag_est(I) = mags;
end