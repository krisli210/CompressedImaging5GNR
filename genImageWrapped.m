function [nse, nse_SS] = genImageWrapped(prm, H_TX, H_RX, Psi_AZ, Psi_R, Phi_R)
    [H_tens, RangeAzProfile, ScatPosPol, threshold] = genGridChannel(prm);
    
    [txGrid, F] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K


    % Rx Signal
    
    % W = exp(-1j*2*pi* ...
    %     sind(randi(359, ([prm.NumRxElements prm.NumRxElements prm.K]))));
    W = repmat(eye(prm.NumRxElements), [1 1 prm.K]);
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    
    for k = 1:prm.K
        Y_tens(:, :, k) = W(:, :, k) * awgn(H_tens(:, :, k) * txGrid(:, :, k), prm.SNR_dB, 'measured');
    end
    
    % % % Receive Processing
    Y_kron = reshape(Y_tens, [prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K]);
    z_theta_per_K = zeros(prm.N_theta, prm.K);

    azSupport = zeros(1, prm.N_theta);

    azSupport_SS = zeros(1, prm.N_theta);
    freqSamples = 1: (floor(prm.K / prm.N_R)) : prm.K;
    numFreqSamples = length(freqSamples);
    % Az Cutting
    for k = 1:prm.K
        Phi_AZ = kron(squeeze(txGrid(:, :, k)).', W(:, :, k));
        [z_theta_per_K(:, k), I] = ...,
        solveCS_OMP(Y_kron(:, k), Phi_AZ, Psi_AZ, 10);
        azSupport(I) = azSupport(I) | 1;

        if(any(k == freqSamples))
            azSupport_SS(I) = azSupport_SS(I) | 1;
        end
    end

    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    RangeAzProfile_hat_SS = zeros(size(RangeAzProfile));
    Phi_R_SS = eye(numFreqSamples);
    for azBin = find(azSupport)
        [z_hat_R, I_R] = solveCS_OMP(z_theta_per_K(azBin, :).', Phi_R, Psi_R);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
        if (azSupport_SS(azBin))
            [z_hat_R_SS, I_R_SS] = solveCS_OMP(z_theta_per_K(azBin, freqSamples).', Phi_R_SS, Psi_R(freqSamples, :));
            RangeAzProfile_hat_SS(:, azBin) = z_hat_R;
        end
    end
    nse = norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2;
    nse_SS = norm(RangeAzProfile_hat_SS - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2;
end