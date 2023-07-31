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
    prm.N_T = 1; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.N_s = prm.N_T * prm.Nofdm;
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
    prm.N_theta = 256;
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    prm.AzBins_Sq = thetaMin:(120/sqrt(prm.N_theta)-1):thetaMax;

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
    %Remove the norm 
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
    prm.Pt_dBm = 20;
    prm.Pt_W = 10^(prm.Pt_dBm/10)*1e-3; % Watts
    
    nIters = 1000;
    SNR_dB_range = -30:5:30;
    SNR_lin_range = 10.^(SNR_dB_range./10);
    container = zeros(3, length(SNR_dB_range), nIters);
    for SNR_ind = 1:length(SNR_dB_range)
        prm.SNR_dB = SNR_dB_range(SNR_ind);
        prm.SNR_lin = SNR_lin_range(SNR_ind);
        for i = 1:nIters
            [container(1, SNR_ind, i), ...
             container(2, SNR_ind, i), ...
             container(3, SNR_ind, i)] = compareW_Wrapped(prm, H_TX, H_RX, Psi_AZ, Psi_R);
        end
    end
    save('compareW_1000_iters', 'container', 'SNR_dB_range');