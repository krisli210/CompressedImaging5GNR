close all
clear 

rng(52);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 30; % number of resource blocks
    prm.K = 12*prm.NRB;
    
    prm.NumUsers = 3; % U per RB
    prm.N_T = 3; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 16; %modulation order
    
% % % END OFDM Signal Params

% % % % % Array and Channel Info

    prm.BsPos = [0; 0; 0];
    prm.BsArraySize = 16; %BS Dimension
    prm.NumBsElements = prod(prm.BsArraySize);
    prm.DeltaTX = .5; % Element spacing normalized by wavelength
    prm.BsAZlim = [-60 60];
    prm.BsELlim = [-90 0];
    
    prm.RxPos = [0; 0; 0];
    prm.RxArraySize = 16;
    prm.DeltaRX = prm.NumBsElements * .5; % Set for virtual array 
    prm.NumRxElements = prod(prm.RxArraySize);
    prm.RxAZlim = prm.BsAZlim;
    prm.RxELlim = [-90 0];
    
    % % % % % % % Target Construction
%     prm.N_theta = prm.BsArraySize * prm.RxArraySize; %number of grid spots, i.e. dimension of the quantized azimuth profile
    prm.N_theta = 128;
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    
    % % % % % % % Grid/Target Construction
    prm.L = 10;
    prm.rMin = 20; prm.rMax = 70;
    
    max_FSPL_dB = 10*log10((4*pi*prm.rMax/prm.lam)^-2);
    limit_rangeRes = false;
    if (limit_rangeRes)
        prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*prm.K); % nominal range resolution
    else
        % Block for changing rangeRes to a non-nominal value
        prm.delta_R = 2;
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
% % % Transmit Signal Construction
    
    %Generates baseband frequency-domain signal across freq-time-space
    [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K

% % % END Transmit Signal Construction
    
    % Rx Signal    
    W = repmat(eye(prm.NumRxElements), [1 1 prm.K]);
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    
    for k = 1:prm.K
        n_k = 0;
        Y_tens(:, :, k) = W(:, :, k) * (H_tens(:, :, k) * txGrid(:, :, k) + n_k);
        % [Y_tens(:, :, k), var] = awgn(H_tens(:, :, k) * txGrid(:, :, k), 100, 'measured') ;
        % Y_tens(:, :, k) = W(:, :, k) * Y_tens(:, :, k);
    end

% % % Receive Processing
    Y_kron = reshape(Y_tens, [prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K]);
    z_theta_per_K = zeros(prm.N_theta, prm.K);

    azSupport = zeros(1, prm.N_theta);

    % Az Cutting
    freqSamples = 1: (floor(prm.K / prm.N_R)) : prm.K;
    % freqSamples = 1:prm.K;
    tic
    for k = freqSamples
        Phi_AZ = kron(squeeze(txGrid(:, :, k)).', W(:, :, k));
        A_Theta = Phi_AZ*Psi_AZ;
        z_theta_per_K(:, k) = omp(A_Theta, Y_kron(:, k), prm.L, 1e-20);
        I = find(z_theta_per_K(:, k));
        azSupport(I) = azSupport(I) | 1;
    end
    toc
    figure; 
    hold on; 
    stem(-60:120/(prm.N_theta-1):60, abs(sum(RangeAzProfile, 1)));
    stem(-60:120/(prm.N_theta-1):60, mean(abs(z_theta_per_K(:, freqSamples)), 2), '--x');
    xlabel('\theta');
    ylabel('$|\hat{\mathbf{z}}_{\Theta}|$', 'Interpreter','latex');
    title('Noisy Azimuth Profile')
    legend({'True', 'Estimate'}, 'Location', 'best')  
    
    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    tic
    for azBin = find(azSupport)
        z_hat_R = omp(Psi_R(freqSamples, :), z_theta_per_K(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    NSE = 10*log10(norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2);
    toc
    figure;

    subplot(1, 2, 1);
    [h, c, lim] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile).^2).', ...
        'typerose', 'default', 'labelR', 'r [m]');
    c.Label.String = 'Measured RCS [dB]';
    sgtitle({'True (L) vs. Estimate (R)', ...
    ['NSE = ' ' $\frac{||\hat{z} - z||^2}{||z||^2} = $' num2str(NSE) ' [dB]']}, ...
    'interpreter', 'latex'); 
    
    subplot(1, 2, 2);
    [h, c_hat, lim] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile_hat).^2).', ...
        'typerose', 'default', 'labelR', 'r [m]', 'lim', lim);
    c_hat.Label.String = 'Measured RCS [dB]';
