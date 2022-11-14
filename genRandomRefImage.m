function [RefImg, Gamma, H_TX, H_RX, physH] = genRandomRefImage(prm, ScatPosPol, ScatPosCart, ScatCoeffs, rMax, thetaMin, thetaMax)
    
    
    if (floor(sqrt(prm.K)) ~= sqrt(prm.K))
        error('Number of voxels should be power of 2')
    end
    % Discretize r\in(rMin, rMax), \theta\in(thetaMin, thetaMax) w/ K voxels
    % Output a sqrt(K) x sqrt(K) image 
    % H_TX returns M x K, H_RX returns N x K 
    
    nScat = size(ScatPosPol, 2);

    rangeBins = 0:rMax/(sqrt(prm.K)):rMax;
    rangeBins = rangeBins(2:end); % JANK!!!
    thetaBins = thetaMin:(thetaMax-thetaMin)/(sqrt(prm.K)-1):thetaMax;
    
    % fill in corresponding bins
    RefImg = zeros(sqrt(prm.K), sqrt(prm.K));
    for scat = 1:nScat
        [~, rangeInd] = min(abs(ScatPosPol(1, scat) - rangeBins));
        [~, thetaInd] = min(abs(ScatPosPol(2, scat) - thetaBins));
        RefImg(rangeInd, thetaInd) = RefImg(rangeInd, thetaInd) + ScatCoeffs(scat); 
    end
    Gamma = reshape(RefImg, [prm.K, 1]);
    diagGamma = diag(Gamma);
    
    rangeOverK = repmat(rangeBins, [1 sqrt(prm.K)]);
    
    H_TX = exp(-1j*(2*pi/prm.lam).*(rangeOverK))./(4*pi.*(rangeOverK)).^2; %Make this variable over BS pos
    H_TX = repmat(H_TX, [prm.NumBsElements 1]);

    H_RX = exp(-1j*(2*pi/prm.lam).*(rangeOverK))./(4*pi.*(rangeOverK)).^2;
    H_RX = repmat(H_RX, [prm.NumRxElements 1]);
    
    physH = H_RX * diagGamma * H_TX.';
end