function [nse] = genImageWrapped(prm, H_TX, H_RX, Psi_AZ, Psi_R, Phi_R)
    [H_tens, RangeAzProfile, ScatPosPol, threshold] = genGridChannel(prm);
    
    [txGrid, F] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K


    % Rx Signal
    SNR_dBm_per_K = prm.Pt + prm.Pr - (prm.No + 10*log10(prm.Delta_f)); % Pt / N_0
    
    % W = exp(-1j*2*pi* ...
    %     sind(randi(359, ([prm.NumRxElements prm.NumRxElements prm.K]))));
    W = repmat(eye(prm.NumRxElements), [1 1 prm.K]);
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    
    for k = 1:prm.K
        n_k = sqrt(10^(-SNR_dBm_per_K / 10)) .* ... 
            (randn([prm.NumRxElements, prm.N_T * prm.Nofdm]) + 1j * randn([prm.NumRxElements, prm.N_T * prm.Nofdm]));
        % n_k = 0;
        Y_tens(:, :, k) = W(:, :, k) * (H_tens(:, :, k) * txGrid(:, :, k) + n_k);
    end
    
    % % % Receive Processing
    Y_kron = reshape(Y_tens, [prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K]);
    z_theta_per_K = zeros(prm.N_theta, prm.K);

    azSupport = zeros(1, prm.N_theta);

    % Az Cutting
    for k = 1:prm.K
        Phi_AZ = kron(squeeze(txGrid(:, :, k)).', W(:, :, k));
        [z_theta_per_K(:, k), I] = ...,
        solveCS_OMP(Y_kron(:, k), Phi_AZ, Psi_AZ, 10);
        azSupport(I) = azSupport(I) | 1;
    end

    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        [z_hat_R, I_R] = solveCS_OMP(z_theta_per_K(azBin, :).', Phi_R, Psi_R);
        % z_hat_R = linsolve(Psi_R, z_theta_per_K(azBin, :).');
        % RangeAzProfile_hat(z_hat_R >= threshold, azBin) = z_hat_R(z_hat_R >= threshold);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    nse = norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2;
end