function [ScatPosPol, ScatPosCart] = generateScatPositions(nScat, rMin, rMax, thetaMin, thetaMax)
    r = (rMax-rMin).*rand(1, nScat) + rMin;
    theta = (thetaMax-thetaMin).*rand(1, nScat) + thetaMin;
    
    ScatPosPol = [r; theta; zeros(1, nScat)];
    ScatPosCart = [r.*cosd(theta) ; r.*sind(theta) ; zeros(1, nScat)];
end