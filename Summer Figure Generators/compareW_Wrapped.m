function [NSE_Y, NSE_MF, NSE_MMSE] = compareW_Wrapped(prm, H_TX, H_RX, Psi_AZ, Psi_R)
    
    [H_tens, RangeAzProfile, ScatPosPol, threshold, a] = genGridChannel(prm);

    [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX, prm.Pt_W); % (Nofdm * N_T) x M x K
    

    freqSamples = 1: (floor(prm.K / prm.N_R)) : prm.K;
    
    z_theta_per_K_wholeY = zeros(prm.N_theta, prm.K);
    z_theta_per_K_MF = zeros(prm.N_theta, prm.K);
    z_theta_per_K_MMSE = zeros(prm.N_theta, prm.K);
    
    azSupport_wholeY = zeros(1, prm.N_theta);
    azSupport_MF = zeros(1, prm.N_theta);
    azSupport_MMSE = zeros(1, prm.N_theta);
    for k = freqSamples
        X_k = txGrid(:, :, k);
        
        % SNR = (Pt / a_k^2) / sigma^2 -> sigma^2 = 
        Y_k = zeros(prm.NumRxElements, prm.N_s);

        P_r = prm.Pt_W / a(k)^2;
        sigma_sq = P_r / (prm.SNR_lin * prm.NumRxElements);
        for n_s = 1:prm.N_s
            n_k = sqrt(sigma_sq /2) * (randn(prm.NumRxElements, 1) + 1j*randn(prm.NumRxElements, 1));
            Y_k(:, n_s) = H_tens(:, :, k) * X_k(:, n_s) + n_k;
        end
        
        W_MF = X_k';
        W_MMSE = pinv(X_k);
        
        H_hat_MF = Y_k * W_MF;
        H_hat_MMSE = Y_k * W_MMSE;

        Y_kron_k = reshape(Y_k, [prm.NumRxElements * prm.N_T * prm.Nofdm, 1]);
        
%         Psi_AZ_normalized = a(k) * Psi_AZ;
        Psi_AZ_normalized = Psi_AZ;

        A_Theta_k = kron(X_k.', eye(prm.NumRxElements)) * Psi_AZ_normalized;
        
        z_theta_per_K_wholeY(:, k) = omp( A_Theta_k, Y_kron_k, prm.L, 1e-20);
        I = find(z_theta_per_K_wholeY(:, k));
        azSupport_wholeY(I) = azSupport_wholeY(I) | 1;

        z_theta_per_K_MF(:, k) = omp( Psi_AZ_normalized, ...
                                 reshape(H_hat_MF, [numel(H_hat_MF), 1]), prm.L, 1e-20);
        I = find(z_theta_per_K_MF(:, k));
        azSupport_MF(I) = azSupport_MF(I) | 1;

        z_theta_per_K_MMSE(:, k) = omp( Psi_AZ_normalized, ...
                                 reshape(H_hat_MMSE, [numel(H_hat_MMSE), 1]), prm.L, 1e-20);
        I = find(z_theta_per_K_MMSE(:, k));
        azSupport_MMSE(I) = azSupport_MMSE(I) | 1;
    
    end
    RangeAzProfile_hat_wholeY = zeros(size(RangeAzProfile));
    for azBin = find(azSupport_wholeY)
        z_hat_R = omp( Psi_R(freqSamples, :), z_theta_per_K_wholeY(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat_wholeY(:, azBin) = z_hat_R;
    end
    NSE_Y = norm(RangeAzProfile_hat_wholeY - RangeAzProfile, 'fro').^2 ./ norm(RangeAzProfile, 'fro').^2;

    RangeAzProfile_hat_MF = zeros(size(RangeAzProfile));
    for azBin = find(azSupport_MF)
        z_hat_R = omp( Psi_R(freqSamples, :), z_theta_per_K_MF(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat_MF(:, azBin) = z_hat_R;
    end
    NSE_MF = norm(RangeAzProfile_hat_MF - RangeAzProfile, 'fro').^2 ./ norm(RangeAzProfile, 'fro').^2;
    
    RangeAzProfile_hat_MMSE = zeros(size(RangeAzProfile));
    for azBin = find(azSupport_MMSE)
        z_hat_R = omp( Psi_R(freqSamples, :), z_theta_per_K_MMSE(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat_MMSE(:, azBin) = z_hat_R;
    end
    NSE_MMSE = norm(RangeAzProfile_hat_MMSE - RangeAzProfile, 'fro').^2 ./ norm(RangeAzProfile, 'fro').^2;

end