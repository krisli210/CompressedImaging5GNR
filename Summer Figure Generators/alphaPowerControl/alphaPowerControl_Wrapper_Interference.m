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
    
    prm.NumUsers = 4; % U per RB
    prm.NumVirtualUsers = 4; % V per RB;
    prm.alpha = 1; % Power scaling. 
    prm.commsSNR_dB = 10;

    prm.N_T = 1; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.N_s = prm.N_T*prm.Nofdm;
    prm.MCS = 16; %modulation order
    
% % % END OFDM Signal Params

% % % % % Array and Channel Inf

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
    prm.L = 133;
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
    for j = 1:prm.N_theta
        H_TX(:, j) = exp(-1j * 2 * pi * prm.DeltaTX * (0:prm.NumBsElements-1) * sind(prm.AzBins(j))).';
        H_RX(:, j) = exp(-1j * 2 * pi * prm.DeltaRX * (0:prm.NumRxElements-1) * sind(prm.AzBins(j))).';
    end
    Psi_AZ = kr(H_TX, H_RX); % 

    Psi_R = zeros(prm.K, prm.N_R);
    for r = 1:length(prm.RangeBins)
        tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
        Psi_R(:, r) = exp(-1j * 2*pi .* [1:prm.K]*prm.Delta_f * tau_r); % This is just the ifft of an identity
    end

    prm.Pt_dBm = 20;
    prm.Pt_W = 10^(prm.Pt_dBm/10)*1e-3; % Watts
    prm.SNR_dB = 20;
    prm.SNR_lin = 10^(prm.SNR_dB/10);
    prm.gamma_dB = 0; %channel gain
    prm.gamma = 10^(prm.gamma_dB/10);
    prm.sigma_N_sq = 1e-6;
    
    nItersPerScene = 50;
    nScenes = 1;
    alpha_range = 1:-.05:0;
    u_range = 0:8;
    v_range = 8 - u_range;
    psnrs = zeros(length(v_range), length(alpha_range), nScenes * nItersPerScene);
    sum_capacities = zeros(length(v_range), length(alpha_range), nScenes * nItersPerScene);
    
    for i = 1:nScenes
        for alpha_ind = 1:length(alpha_range)
            prm.alpha = alpha_range(alpha_ind);
            for v_ind = 1:length(v_range)
                prm.NumVirtualUsers = v_range(v_ind);
                prm.NumUsers = u_range(v_ind);
                [H_tens, RangeAzProfile, ~, ~, a] = genGridChannel(prm);
                for j = 1:nItersPerScene
                    currVecInd = (i-1)*nItersPerScene + j;
                    [sum_capacities(v_ind, alpha_ind, currVecInd), psnrs(v_ind, alpha_ind, currVecInd)] = alphaPowerControlWrapped_Interference(prm, H_TX, H_RX, Psi_AZ, Psi_R, H_tens, RangeAzProfile, a);
                end
            end
        end
    end

    save("C:\Users\krisl\Desktop\Summer2023\CompressedImaging5GNR\Summer Figure Generators\alphaPowerControl\alphaPowerControl_vary_V_RICEscene_Interference_v3_0dBGain_fullData", ...
        'psnrs', 'sum_capacities', 'prm', 'alpha_range', 'v_range', 'u_range');