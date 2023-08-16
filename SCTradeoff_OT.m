close all
clear 

rng(42);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 60; % number of resource blocks
    prm.K = 12*prm.NRB;
    
    prm.NumUsers = 3; % U per RB
    prm.NumVirtualUsers = 5; % V per RB;
    prm.alpha = 1; % Power scaling. 
    prm.N_T = 1; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.N_s = prm.N_T*prm.Nofdm;
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
    prm.L = 30;
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

    [H_tens, RangeAzProfile, ScatPosPol, threshold, a] = genGridChannel(prm);
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
    % Generate the user location distribution
    [theta_dist, theta_dist_cum, theta_dist_v, theta_dist_cum_v] = getUserDistribution(prm.N_theta, prm.AzBins, 'custom bimodal');
    
    %Generates baseband frequency-domain signal across freq-time-space
    prm.Pt_dBm = 20;
    prm.Pt_W = 10^(prm.Pt_dBm/10)*1e-3; % Watts
    % [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX, prm.Pt_W); % (Nofdm * N_T) x M x K
    [txGrid, powerLog, userAngleIndLog] = genAngleDefFreqTxGrid_v2(prm.NumBsElements, prm.NumUsers, prm.NumVirtualUsers, prm.alpha, prm.MCS, ...
                prm.N_T, prm.Nofdm, prm.K, H_TX, prm.Pt_W, theta_dist_cum, theta_dist_cum_v);
% % % END Transmit Signal Construction
    
    % Rx Signal    
    prm.SNR_dB = 20;
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
        
        z_theta_per_K(:, k) = omp( A_Theta_k, Y_kron_k, prm.L, 1e-20);
        I = find(z_theta_per_K(:, k));
        azSupport(I) = azSupport(I) | 1;
    end
    % figure; 
    % hold on; 
    % stem(-60:120/(prm.N_theta-1):60, abs(sum(RangeAzProfile, 1)));
    % stem(-60:120/(prm.N_theta-1):60, mean(abs(z_theta_per_K(:, 1)), 2), '--x');
    % xlabel('\theta');
    % ylabel('$|\hat{\mathbf{z}}_{\Theta}|$', 'Interpreter','latex');
    % title('Noisy Azimuth Profile')
    % legend({'True', 'Estimate'}, 'Location', 'best')  
    
    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        z_hat_R = omp(Psi_R(freqSamples, :), z_theta_per_K(azBin, freqSamples).', 10, 1e-20);
        RangeAzProfile_hat(:, azBin) = z_hat_R;
    end
    NSE = 10*log10(norm(RangeAzProfile_hat - RangeAzProfile).^2 ./ norm(RangeAzProfile).^2);
    peaksnr = psnr(abs(RangeAzProfile_hat), abs(RangeAzProfile), max(abs(RangeAzProfile), [], "all"));

    % figure;
    % subplot(1, 2, 1);
    % [h, c, lim] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile).^2).', ...
    %     'typerose', 'default', 'labelR', 'r [m]');
    % c.Label.String = 'Measured RCS [dB]';
    % sgtitle({'True (L) vs. Estimate (R)', ...
    % ['NSE = ' ' $\frac{||\hat{z} - z||^2}{||z||^2} = $' num2str(NSE) ' [dB]']}, ...
    % 'interpreter', 'latex'); 
    
    % subplot(1, 2, 2);
    % [h, c_hat, lim] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile_hat).^2).', ...
    %     'typerose', 'default', 'labelR', 'r [m]', 'lim', lim);
    % c_hat.Label.String = 'Measured RCS [dB]';
    
    figure; hold on;
    plot(prm.AzBins, theta_dist);
    plot(prm.AzBins, theta_dist_v);
    xlabel('Azimuth [$^\circ$]', 'Interpreter','latex', 'FontSize',14);
    ylabel('Probability of User Angle')
    title('PMF of User Locations')
    legend({'Real', 'Virtual'}, 'Location','best', 'FontSize', 14)

    figure;
    imagesc(prm.AzBins, prm.RangeBins, abs(RangeAzProfile));
    xlabel('Azimuth [$^\circ$]', 'Interpreter','latex', 'FontSize',14);
    ylabel('Range $[m]$', 'Interpreter','latex', 'FontSize',14)
    title('True Range-Angle Profile')

    figure;
    imagesc(prm.AzBins, prm.RangeBins, abs(RangeAzProfile_hat));
    xlabel('Azimuth [$^\circ$]', 'Interpreter','latex', 'FontSize',14);
    ylabel('Range $[m]$', 'Interpreter','latex', 'FontSize',14)
    title(['Image Reconstruction, PSNR = ', num2str(peaksnr), ' [dB]'], ...
           'Interpreter','latex','FontSize',14)

    % sgtitle({'True (L) vs. Estimate (R)', ...
    % ['NSE = ' num2str(NSE) ' [dB]']}, ...
    % 'interpreter', 'latex'); 
