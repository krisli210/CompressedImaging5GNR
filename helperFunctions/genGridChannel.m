function [H_tens, RangeAzProfile, ScatPosPol, threshold] = genGridChannel(prm)
    % H _tens output as N x M x K
    % RangeAzProfile as N_R x N_theta
    % ScatPosPol as 3 x L

%     azInd = randperm(prm.N_theta, prm.L);
%     rangeInd = randperm(prm.N_R, prm.L);
    % azInd = randi(prm.N_theta, [prm.L, 1]);
    % rangeInd = randi(prm.N_R, [prm.L, 1]);
    [azInd, rangeInd] = getRICE_Image();
    prm.L = length(azInd);

    azValues = prm.AzBins(azInd);
    rangeValues = prm.RangeBins(rangeInd);
    ScatPosPol = [rangeValues; azValues; zeros(1, prm.L)];
    ScatCoeff = ones(1, prm.L) .* complex(1, 1) ./ sqrt(2); %Unit reflectors

    RangeAzProfile = zeros(prm.N_R, prm.N_theta);
    H_tens = zeros(prm.NumRxElements, prm.NumBsElements, prm.K);
    
    for l = 1:prm.L
        tau_r = 2*rangeValues(l) / prm.PropagationSpeed;
        tau_m = prm.DeltaTX * (0:prm.NumBsElements-1)*sind(azValues(l)); %Check this
        tau_n = prm.DeltaRX * (0:prm.NumRxElements-1)*sind(azValues(l));
        tau_n_m = zeros(prm.NumRxElements, prm.NumBsElements);

        for n = 1:prm.NumRxElements
            for m = 1:prm.NumBsElements
                tau_n_m(n,m) = tau_n(n) + tau_m(m);
            end
        end

        PL = (4*pi * 2*rangeValues(l)/prm.lam)^-2;
        RangeAzProfile(rangeInd(l), azInd(l)) = PL*ScatCoeff(l);
        for k = 1:prm.K
            % Currently at the actual frequency content of k rather than
            % just k-th bin - this varies between models I've seen

            %In Guan, the baseband received signal has the k-th subcarrier
            %scaled by the SCS, while in Araujo (8) it is only at the bin
            H_tens(:, :, k) = H_tens(:, :, k) + PL*ScatCoeff(l)*exp(-1j * 2*pi * ( (k*prm.Delta_f*tau_r) + tau_n_m) ); % Is the az-induced delay scaled by freq? - NO
        end
    end

    threshold = complex(1, 1) ./ sqrt(2) * (4*pi* 2 * (prm.rMax+10)/prm.lam)^-2; % Threshold is given as a unit reflector past max range
end