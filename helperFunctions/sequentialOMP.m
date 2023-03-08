function RangeAzProfileHat = sequentialOMP(txGrid, Y_tens, Psi_AZ, Psi_R, H_TX, H_RX, prm)
    nIters = prm.L; % for now
    RangeAzProfileHat = zeros(prm.N_R, prm.N_theta);

    W = eye(prm.NumRxElements);
    for i = 1:nIters
        freqSample = randperm(prm.K, 1); 
        x_AZ = txGrid(:, :, 1); % Space x Time - N x Nofdm * N_T
        y_AZ = Y_tens(:, :, 1);

        y_vec_AZ = reshape(y_AZ, [numel(y_AZ), 1]);
        Phi_AZ = kron(x_AZ.', W); % N*Nofdm*N_T x MN

        [z_hat_AZ, I_AZ] = solveCS_OMP(y_vec_AZ, Phi_AZ, Psi_AZ, 1);

%         tx_rand = randi(prm.BsArraySize);
%         rx_rand = randi(prm.RxArraySize);
    
%         x_R = squeeze(txGrid(tx_rand, :, :)).'; % Freq x Time
%         y_R = squeeze(Y_tens(rx_rand, :, :)).'; % 
        weight_mags = ones(1, length(I_AZ)).';
        y_R = steerTensor(Y_tens, H_RX, I_AZ, weight_mags).';
        x_R = steerTensor(txGrid, H_TX, I_AZ, weight_mags).';
        
        Phi_R = eye(size(Psi_R, 1));
        H_hat = mean(y_R ./ x_R, 2); % Average freq response over OFDM symbols
        [z_hat_R, I_R] = solveCS_OMP(H_hat, Phi_R, Psi_R, 1);
        
        RangeAzProfileHat(I_R, I_AZ) = z_hat_AZ(I_AZ);
%         RangeAzProfileHat(:, I_AZ) = z_hat_R;
        Psi_AZ(:, I_AZ) = 0;
    end
end