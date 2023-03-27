close all
clear 

rng(42);


% % % OFDM Signal Params

    prm.CenterFreq = 26e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 30; % number of resource blocks
    prm.K = 12*prm.NRB;
    
    prm.NumUsers = 4; % U per RB
    prm.N_T = 1; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 16; %modulation order
    
    Pt = 20; %dBm 
    Pr = 20; %dBm

    No = -174; %dBm/Hz
% % % END OFDM Signal Params

% % % % % Array and Channel Info

    prm.BsPos = [0; 0; 0];
    prm.BsArraySize = 18; %BS Dimension
    prm.NumBsElements = prod(prm.BsArraySize);
    prm.DeltaTX = .5; % Element spacing normalized by wavelength
    prm.BsAZlim = [-60 60];
    prm.BsELlim = [-90 0];
    
    prm.RxPos = [0; 0; 0];
    prm.RxArraySize = 15;
    prm.DeltaRX = prm.NumBsElements * .5; % Set for virtual array 
    prm.NumRxElements = prod(prm.RxArraySize);
    prm.RxAZlim = prm.BsAZlim;
    prm.RxELlim = [-90 0];
    
    % % % % % % % Target Construction
%     prm.N_theta = prm.BsArraySize * prm.RxArraySize; %number of grid spots, i.e. dimension of the quantized azimuth profile
    prm.N_theta = 24;
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    
    % % % % % % % Grid/Target Construction
    prm.L = 5;
    prm.rMin = 20; prm.rMax = 70;
    
    max_FSPL_dB = 10*log10((4*pi*prm.rMax/prm.lam)^-2);
    limit_rangeRes = true;
    if (limit_rangeRes)
        prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*prm.K); % nominal range resolution
    else
        % Block for changing rangeRes to a non-nominal value
        prm.delta_R = 5;
    end
    prm.WholeRange = 0:prm.delta_R:(prm.K-1)*prm.delta_R;
    minIndex = find(prm.WholeRange < prm.rMin, 1, 'last')+1;
    maxIndex = find(prm.WholeRange > prm.rMax, 1, 'first')-1;

%     prm.RangeBins = prm.rMin : prm.delta_R : prm.rMax;
    prm.RangeBins = prm.WholeRange(minIndex:maxIndex);
    prm.N_R = numel(prm.RangeBins);

    [H_tens, RangeAzProfile, ScatPosPol, threshold] = genGridChannel(prm);
    % % % % % % % END Grid/Target Construction

    % Spatial Dictionary Construction
    H_TX = zeros(prm.NumBsElements, prm.N_theta);
    H_RX = zeros(prm.NumRxElements, prm.N_theta);
    for n = 1:prm.N_theta
        H_TX(:, n) = exp(-1j * 2 * pi * prm.DeltaTX * (0:prm.NumBsElements-1) * sind(prm.AzBins(n))).';
        H_RX(:, n) = exp(-1j * 2 * pi * prm.DeltaRX * (0:prm.NumRxElements-1) * sind(prm.AzBins(n))).';
    end
    Psi_AZ = kr(H_TX, H_RX); % 

    Psi_R = zeros(prm.K, prm.N_R);
    for r = 1:length(prm.RangeBins)
        tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
        Psi_R(:, r) = exp(-1j * 2*pi .* [1:prm.K]*prm.Delta_f * tau_r); % This is just the ifft of an identity
    end
    Phi_R = eye(size(Psi_R, 1));
% % % Transmit Signal Construction
    
    %Generates baseband frequency-domain signal across freq-time-space
    %Both a randomly precoded and ideally precoded/data structure 
    [txGrid, F] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K
    
% % % END Transmit Signal Construction
    
    % Rx Signal
    SNR_dB_per_K = Pt + Pr - (No + 10*log10(prm.Delta_f)); % Pt / N_0
    
    % W = exp(-1j*2*pi* ...
    %     sind(randi(359, ([prm.NumRxElements prm.NumRxElements prm.K]))));
    W = repmat(eye(prm.NumRxElements), [1 1 prm.K]);
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    
    for k = 1:prm.K
        n_k = sqrt(10^(-SNR_dB_per_K / 10)) .* ... 
            (randn([prm.NumRxElements, prm.N_T * prm.Nofdm]) + 1j * randn([prm.NumRxElements, prm.N_T * prm.Nofdm]));
        n_k = 0;
        Y_tens(:, :, k) = W(:, :, k) * (H_tens(:, :, k) * txGrid(:, :, k) + n_k);
    end

% % % Receive Processing
    Y_kron = reshape(Y_tens, [prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K]);
    z_theta_per_K = zeros(prm.N_theta, prm.K);

    azSupport = zeros(1, prm.N_theta);

    % Az Cutting
    for k = 1:prm.K
        Phi_AZ = kron(squeeze(txGrid(:, :, k)).', W(:, :, k));
        [z_theta_per_K(:, k), I] = solveCS_OMP(Y_kron(:, k), Phi_AZ, Psi_AZ, 10);
        % z_theta_per_K(:, k) = linsolve(Phi_AZ*Psi_AZ, Y_kron(:, k));
        % I = find(z_theta_per_K(:, k));
        azSupport(I) = azSupport(I) | 1;
    end

    figure; 
    subplot(1, 2, 1); hold on; 
    stem(-60:120/(prm.N_theta-1):60, abs(sum(RangeAzProfile, 1)));
    % stem(-60:120/(prm.N_theta-1):60, mean(abs(z_theta_per_K), 2), '--x');
    stem(-60:120/(prm.N_theta-1):60, abs(z_theta_per_K(:, 1)), '--x');
    xlabel('\theta');
    ylabel('Magnitude');
    title('Azimuth Profile')
    legend({'True', 'Estimate'}, 'Location', 'south')  
    
    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        [z_hat_R, I_R] = solveCS_OMP(z_theta_per_K(azBin, :).', Phi_R, Psi_R);
        % z_hat_R = linsolve(Psi_R, z_theta_per_K(azBin, :).');
        % RangeAzProfile_hat(z_hat_R >= threshold, azBin) = z_hat_R(z_hat_R >= threshold);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    ssim_val = ssim(abs(RangeAzProfile_hat).^2, ...
        abs(RangeAzProfile).^2);
    NSE = 10*log10(norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2);

    figure;
    subplot(1, 2, 2);
    [h, c_hat, lim] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile_hat).^2).', ...
        'typerose', 'default', 'labelR', 'r [m]');
    c_hat.Label.String = 'Measured RCS [dB]';
    subplot(1, 2, 1);

    [h, c, lim] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile).^2).', ...
        'typerose', 'default', 'labelR', 'r [m]', 'lim', lim);
    c.Label.String = 'Measured RCS [dB]';
    sgtitle({'True (L) vs. Estimate (R)', ...
    ['NSE = ' ' $\frac{||\hat{z} - z||^2}{||z||^2} = $' num2str(NSE) ' [dB]']}, ...
    'interpreter', 'latex'); 
