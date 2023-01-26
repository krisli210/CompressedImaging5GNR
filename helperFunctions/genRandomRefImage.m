function [RefImg, Gamma, H_TX, H_RX, physH] = genRandomRefImage(prm, ScatPosCart, ScatCoeffs, BsArrayPos, RxArrayPos, rMax, thetaMin, thetaMax)
    
    
    if (floor(sqrt(prm.K)) ~= sqrt(prm.K))
        error('Number of voxels should be power of 2')
    end
    % Discretize r\in(rMin, rMax), \theta\in(thetaMin, thetaMax) w/ K voxels
    % Output a sqrt(K) x sqrt(K) image 
    % H_TX returns M x K, H_RX returns N x K 
    
    nScat = size(ScatPosCart, 2);

    rangeBins = 0:rMax/(sqrt(prm.K)):rMax;
    rangeBins = rangeBins(2:end); % JANK!!!
    thetaBins = thetaMin:(thetaMax-thetaMin)/(sqrt(prm.K)-1):thetaMax;
    
    RefImgCoorPol = [repmat(rangeBins, [1 sqrt(prm.K)]); reshape(repmat(thetaBins, [sqrt(prm.K) 1]), [1 prm.K]); zeros(1, prm.K)];

    [x, y, z] = pol2cart(deg2rad(RefImgCoorPol(2, :)), RefImgCoorPol(1, :), RefImgCoorPol(3, :));
    RefImgCoorCart = [x; y; z];

    % fill in corresponding bins
    Gamma = zeros(prm.K, 1);
    for scat = 1:nScat
        [~, bin] = min(vecnorm(ScatPosCart(:, scat) - RefImgCoorCart, 2, 1));
        Gamma(bin) = Gamma(bin) + ScatCoeffs(scat);
    end
    
    RefImg = reshape(Gamma, [sqrt(prm.K), sqrt(prm.K)]);
    
    H_TX = zeros(prm.NumBsElements, prm.K);
    H_RX = zeros(prm.NumRxElements, prm.K);
    
    for k = 1:prm.K
        d_TX_K = vecnorm(BsArrayPos - RefImgCoorCart(:, k), 2, 1);
        d_K_RX = vecnorm(RxArrayPos - RefImgCoorCart(:, k), 2, 1);

        H_TX(:, k) = exp(-1j*(2*pi/prm.lam).*(d_TX_K))./(4*pi.*(d_TX_K)).^2;
        H_RX(:, k) = exp(-1j*(2*pi/prm.lam).*(d_K_RX))./(4*pi.*(d_K_RX)).^2;
    end
    
%     %Randomly 0 receiver elements as a test to reduce coherence
%     %i.e., 0 along the rows of H_RX
%     p = .2; % probability we use
%     rows = sort(randperm(prm.NumRxElements, floor(prm.NumRxElements * p)));
%     H_RX = H_RX(rows, :);
%  
%     % DID NOT WORK 

    physH = H_RX * diag(Gamma) * H_TX.';
end