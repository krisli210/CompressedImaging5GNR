function [ScatPosPol, ScatPosCart, azProfile] = generateDiscreteScatPositions(nScat, rMin, rMax, AzBins, ScatCoeffs)
    r = (rMax-rMin).*rand(1, nScat) + rMin;
    
    thetaInd = randperm(length(AzBins), nScat);

    ScatPosPol = [r; AzBins(thetaInd); zeros(1, nScat)];
    ScatPosCart = [r.*cosd(AzBins(thetaInd)) ; r.*sind(AzBins(thetaInd)) ; zeros(1, nScat)];
    
    azProfile = zeros(length(AzBins), 1);
    azProfile(thetaInd) = ScatCoeffs;
end