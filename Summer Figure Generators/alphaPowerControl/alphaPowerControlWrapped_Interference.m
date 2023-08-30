function [sum_capacity, peaksnr] = alphaPowerControlWrapped_Interference(prm, H_TX, H_RX, Psi_AZ, Psi_R, H_tens, RangeAzProfile, a)
    
    % [H_tens, RangeAzProfile, ~, ~, a] = genGridChannel(prm);

    [theta_dist, theta_dist_cum, theta_dist_v, theta_dist_cum_v] = getUserDistribution(prm.N_theta, prm.AzBins, 'custom bimodal');
    theta_dist_cum(end) = 1; %enforcing this cuz rounding errors cause this to be < 1 and destroy the find() on rand
    theta_dist_cum_v(end) = 1;

    %Generates baseband frequency-domain signal across freq-time-space
    prm.Pt_dBm = 20;
    prm.Pt_W = 10^(prm.Pt_dBm/10)*1e-3; % Watts
    % [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX, prm.Pt_W); % (Nofdm * N_T) x M x K
    prm.commsSNR_lin = 10^(prm.commsSNR_dB/10);
        [txGrid, interferenceLog, rateLog] = genAngleDefFreqTxGrid_v2(prm.NumBsElements, prm.NumUsers, prm.NumVirtualUsers, prm.alpha, prm.MCS, ...
                prm.N_T, prm.Nofdm, prm.K, H_TX, prm.Pt_W, theta_dist_cum, theta_dist_cum_v, prm.gamma, prm.sigma_N_sq);

    sum_capacity = mean(sum(rateLog, 1));
% % % END Transmit Signal Construction
    
    % Rx Signal    
    prm.SNR_lin = 10^(prm.SNR_dB/10); 
    W = zeros(prm.NumRxElements, prm.NumRxElements, prm.K);
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    freqSamples = 1: (floor(prm.K / prm.N_R)) : prm.K;

    z_theta_per_K = zeros(prm.N_theta, prm.K);
    azSupport = zeros(1, prm.N_theta);
    for k = freqSamples
        X_k = txGrid(:, :, k);
        
        Y_k = zeros(prm.NumRxElements, prm.N_s);
        
        P_r = prm.Pt_W / a(k)^2; % norms defnied in loyka #2 
        sigma_sq = P_r / (prm.SNR_lin * prm.NumRxElements);
        for n_s = 1:prm.N_s
            n_k = sqrt(sigma_sq /2) * (randn(prm.NumRxElements, 1) + 1j*randn(prm.NumRxElements, 1));
            Y_k(:, n_s) = H_tens(:, :, k) * X_k(:, n_s) + n_k;
        end
        
        Y_kron_k = reshape(Y_k, [prm.NumRxElements * prm.N_T * prm.Nofdm, 1]);
        
        A_Theta_k = kron(X_k.', eye(prm.NumRxElements)) * Psi_AZ;
        
        z_theta_per_K(:, k) = omp( A_Theta_k, Y_kron_k, 30, 1e-20);
        I = find(z_theta_per_K(:, k));
        azSupport(I) = azSupport(I) | 1;
    end
    
    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        z_hat_R = omp(Psi_R(freqSamples, :), z_theta_per_K(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    NSE = 10*log10(norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2);
    peaksnr = psnr(abs(RangeAzProfile_hat), abs(RangeAzProfile), max(abs(RangeAzProfile), [], "all"));

end