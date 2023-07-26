function [nse] = genImageWrapped_SS_randW(prm, H_TX, Psi_AZ, Psi_R, W_setting)
    [H_tens, RangeAzProfile] = genGridChannel(prm);
    
    [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K


    % Rx Signal
    switch W_setting
        case 'eye'
            W = repmat(eye(prm.NumRxElements), [1 1 prm.K]);
        case 'rand'
            W = (2*rand(([prm.NumRxElements prm.NumRxElements prm.K]))-1) .* ...
                exp(-1j*2*pi*sind(randi(359, ([prm.NumRxElements prm.NumRxElements prm.K]))))./sqrt(prm.NumRxElements);
        case 'ffteye'
            W = repmat(fft(eye(prm.NumRxElements)), [1 1 prm.K])./sqrt(prm.NumRxElements);
    end
    
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    
    freqSamples = 1 : (floor(prm.K / prm.N_R)) : prm.K;
    for k = freqSamples
        Y_tens(:, :, k) = W(:, :, k) * awgn(H_tens(:, :, k) * txGrid(:, :, k), prm.SNR_dB, 'measured');
        % Y_tens(:, :, k) = W(:, :, k) * H_tens(:, :, k) * txGrid(:, :, k); 
    end
    
    % % % Receive Processing
    Y_kron = reshape(Y_tens, [prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K]);
    z_theta_per_K = zeros(prm.N_theta, prm.K);

    azSupport = zeros(1, prm.N_theta);

    % Az Cutting
    for k = freqSamples
        Phi_AZ = kron(squeeze(txGrid(:, :, k)).', W(:, :, k));
        z_theta_per_K(:, k) = omp(Phi_AZ*Psi_AZ, Y_kron(:, k), prm.L, 1e-20);
        I = find(z_theta_per_K(:, k));
        azSupport(I) = azSupport(I) | 1;
    end

    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        z_hat_R = omp(Psi_R(freqSamples, :), z_theta_per_K(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    nse = norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2;
end