function [mag_est, I, YI, R] = solveCS_OMP(y, Phi, Psi, varargin)
    if (nargin == 4)
        OMP_iters = varargin{1};
    else
        OMP_iters = 10;
    end
    
    %returns estimated magnitudes mag_est under support I 
    A = Phi * Psi;
    sensingDict = sensingDictionary('CustomDictionary', A);
    [mag_est, YI, I, R] = matchingPursuit(sensingDict, y, maxIterations=OMP_iters, Algorithm="OMP", maxerr={"L1", 1e-4});
    
    A_sub = Phi * Psi(:, I); % Solve magnitude posthence via direct linsolve against estimated support
    if length(I) > 1
        mags = linsolve(A_sub, y);
    else
        mags = linsolve(A_sub, y.');
    end
    mag_est(I) = mags;
end