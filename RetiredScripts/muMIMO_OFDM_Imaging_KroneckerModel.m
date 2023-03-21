close all
clear 

rng(422);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 30; % number of resource blocks
    prm.K = 12*prm.NRB;
    
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
%     prm.DeltaRX = .5;
    prm.NumRxElements = prod(prm.RxArraySize);
    prm.RxAZlim = prm.BsAZlim;
    prm.RxELlim = [-90 0];
    
    %Arrays as uniform linear given in PA toolbox
    BsArray = phased.ULA(prm.NumBsElements, prm.DeltaTX*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true), 'ArrayAxis', 'y');
    
    RxArray = phased.ULA(prm.NumRxElements, prm.DeltaRX*prm.lam, 'Element', phased.IsotropicAntennaElement,'ArrayAxis', 'y');
    
    % % % % % % % Target Construction
%     prm.N_theta = prm.BsArraySize * prm.RxArraySize; %number of grid spots, i.e. dimension of the quantized azimuth profile
    prm.N_theta = 32;
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    
    % % % % % % % Grid/Target Construction
    prm.L = 4;
    prm.rMin = 20; prm.rMax = 100;

    limit_rangeRes = true;
    if (limit_rangeRes)
        prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*prm.K); % nominal range resolution
    else
        % Block for changing rangeRes to a non-nominal value
        prm.delta_R = 6;
    end
    prm.WholeRange = 0:prm.delta_R:(prm.K-1)*prm.delta_R;
    minIndex = find(prm.WholeRange < prm.rMin, 1, 'last')+1;
    maxIndex = find(prm.WholeRange > prm.rMax, 1, 'first')-1;

%     prm.RangeBins = prm.rMin : prm.delta_R : prm.rMax;
    prm.RangeBins = prm.WholeRange(minIndex:maxIndex);
    prm.N_R = numel(prm.RangeBins);

    [H_tens, RangeAzProfile, ScatPosPol, azInd] = genGridChannel(prm);
    % % % % % % % END Grid/Target Construction

    % Spatial Dictionary Construction
    H_TX = zeros(prm.NumBsElements, prm.N_theta);
    H_RX = zeros(prm.NumRxElements, prm.N_theta);
    for n = 1:prm.N_theta
        H_TX(:, n) = (1/sqrt(prm.NumBsElements)) * exp(-1j * 2 * pi * prm.DeltaTX * (0:prm.NumBsElements-1) * sind(prm.AzBins(n))).';
        H_RX(:, n) = (1/sqrt(prm.NumRxElements)) * exp(-1j * 2 * pi * prm.DeltaRX * (0:prm.NumRxElements-1) * sind(prm.AzBins(n))).';
%         H_TX(:, n) = collectPlaneWave(BsArray, 1, [prm.AzBins(n), 0].', prm.CenterFreq);
%         H_RX(:, n) = collectPlaneWave(RxArray, 1, [prm.AzBins(n), 0].', prm.CenterFreq);
    end
    
    Psi_AZ = kr(H_TX, H_RX); % 
    Psi_AZ = Psi_AZ ./ norm(Psi_AZ);

    Psi_R = zeros(prm.K, prm.N_R);
    for r = 1:length(prm.RangeBins)
        tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
        Psi_R(:, r) = exp(-1j * 2*pi .* [0:prm.K-1]*prm.Delta_f * tau_r); % This is just the ifft of an identity
%         Psi_R(:, r) = Psi_R(:, r) ./ norm(Psi_R(:, r));
    end
    Psi_R = Psi_R ./ norm(Psi_R);
% % % Transmit Signal Construction

    prm.NumUsers = 4;
    prm.N_T = prm.NRB; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 16; %modulation order
    
    %Generates baseband frequency-domain signal across freq-time-space
    %Both a randomly precoded and ideally precoded/data structure 
    [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K
    
% % % END Transmit Signal Construction
    
    % Rx Signal
    rxGain = 10^(0/10);
    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    for k = 1:prm.K
%         Y_tens(:, :, k) = rxGain .* (H_tens(:, :, k) * txGrid(:, :, k));
        Y_tens(:, :, k) = rxGain .* (H_tens(:, :, k) * txGrid(:, :, k)); % this kronecker model will only work for now if we use same transmit signal across all freq
    end
    
% % % Receive Processing
    Y_kron = reshape(Y_tens, [prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K]);
    z_theta_per_K = zeros(prm.N_theta, prm.K);
    YI_per_K = zeros(prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K);
    R_per_K = zeros(prm.NumRxElements * prm.N_T * prm.Nofdm, prm.K);

    azSupport = zeros(1, prm.N_theta);
    W = eye(prm.NumRxElements);

    % Az Cutting
    for k = 1:prm.K
        Phi_theta = kron(squeeze(txGrid(:, :, k)).', W);
        [z_theta_per_K(:, k), I, YI_per_K(:, k), R_per_K(:,k)] = ...,
        solveCS_OMP(Y_kron(:, k), Phi_theta, Psi_AZ, 10);
        azSupport(I) = azSupport(I) | 1;
    end

    figure; 
    subplot(1, 2, 1); hold on; 
    stem(-60:120/(prm.N_theta-1):60, abs(sum(RangeAzProfile, 1)));
    stem(-60:120/(prm.N_theta-1):60, abs(mean(z_theta_per_K, 2)), '--x');
    xlabel('\theta');
    ylabel('Magnitude');
    title('Azimuth Profile')
    legend({'True', 'Estimate'}, 'Location', 'south')  
    
    RangeAzProfile_hat = zeros(size(RangeAzProfile));
    for azBin = find(azSupport)
        rangeProfile = ifft(z_theta_per_K(azBin, :));
        RangeAzProfile_hat(:, azBin) = rangeProfile(minIndex:maxIndex);
    end

    figure;
    subplot(1, 2, 1);
    [h, c] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile).^2).', ...
        'typerose', 'default', 'labelr', 'r [m]');
    c.Label.String = 'Measured Reflection Power [dB]';
    
    subplot(1, 2, 2);
    [h, c_hat] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile_hat).^2).', ...
        'typerose', 'default', 'labelr', 'r [m]');
%     c.Label.String = 'Measured Reflection Power [dB]';
    c_hat.Limits = c.Limits;
    sgtitle('True (L) vs. Estimate (R)'); 