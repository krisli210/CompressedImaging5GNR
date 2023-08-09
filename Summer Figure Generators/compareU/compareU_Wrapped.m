function [NSE] = compareU_Wrapped(prm, H_TX, H_RX, Psi_AZ, Psi_R)
    
    [H_tens, RangeAzProfile, ScatPosPol, threshold, a] = genGridChannel(prm);

    [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX, prm.Pt_W); % (Nofdm * N_T) x M x K
    

    freqSamples = 1: (floor(prm.K / prm.N_R)) : prm.K;
    
    z_theta_per_K = zeros(prm.N_theta, prm.K);
    azSupport = zeros(1, prm.N_theta);
    for k = freqSamples
        X_k = txGrid(:, :, k);
        
        % SNR = (Pt / a_k^2) / sigma^2 -> sigma^2 = 
        Y_k = zeros(prm.NumRxElements, prm.N_s);
        
        P_r = prm.Pt_W / a(k)^2; % norms defnied in loyka #2 
        sigma_sq = P_r / (prm.SNR_lin * prm.NumRxElements);
        for n_s = 1:prm.N_s
            n_k = sqrt(sigma_sq /2) * (randn(prm.NumRxElements, 1) + 1j*randn(prm.NumRxElements, 1));
            Y_k(:, n_s) = H_tens(:, :, k) * X_k(:, n_s) + n_k;
        end
        
        Y_kron_k = reshape(Y_k, [prm.NumRxElements * prm.N_T * prm.Nofdm, 1]);
        
%         Psi_AZ_normalized = a(k) * Psi_AZ_unnormalized;
        Psi_AZ_normalized = Psi_AZ;
        A_Theta_k = kron(X_k.', eye(prm.NumRxElements)) * Psi_AZ_normalized;
        
        z_theta_per_K(:, k) = omp( A_Theta_k, Y_kron_k, prm.L, 1e-20);
        I = find(z_theta_per_K(:, k));
        azSupport(I) = azSupport(I) | 1;
    end
    
    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        z_hat_R = omp( Psi_R(freqSamples, :), z_theta_per_K(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    NSE = norm(RangeAzProfile_hat - RangeAzProfile, 'fro').^2 ./ norm(RangeAzProfile, 'fro').^2;

end