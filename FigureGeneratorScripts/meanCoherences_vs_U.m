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
    
    W = repmat(eye(prm.NumRxElements), [1 1 prm.K]);
    % % % % % % % Target Construction
%     prm.N_theta = prm.BsArraySize * prm.RxArraySize; %number of grid spots, i.e. dimension of the quantized azimuth profile
    N_Theta_range = [64:64:320];
    U_range = [1:8];
    
    nIters = 100;
    MCs = zeros(length(N_Theta_range), length(U_range), nIters);
    for N_Theta_ind = 1:length(N_Theta_range)
        prm.N_theta = N_Theta_range(N_Theta_ind);
        thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
        prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    
        % Spatial Dictionary Construction
        H_TX = zeros(prm.NumBsElements, prm.N_theta);
        H_RX = zeros(prm.NumRxElements, prm.N_theta);
        for n = 1:prm.N_theta
            H_TX(:, n) = exp(-1j * 2 * pi * prm.DeltaTX * (0:prm.NumBsElements-1) * sind(prm.AzBins(n))).';
            H_RX(:, n) = exp(-1j * 2 * pi * prm.DeltaRX * (0:prm.NumRxElements-1) * sind(prm.AzBins(n))).';
        end
        Psi_AZ = kr(H_TX, H_RX); % 

        for U_ind = 1:length(U_range)
            prm.NumUsers = U_range(U_ind);
            for n = 1:nIters
                
                % % % Transmit Signal Construction
                    
                    %Generates baseband frequency-domain signal across freq-time-space
                    [txGrid] = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K
                
                % % % END Transmit Signal Construction
                
                % % % Receive Processing
                    
                    azSupport = zeros(1, prm.N_theta);
                    
                    % Az Cutting
                    freqSamples = 1:prm.K;
                    run_MCs = zeros(1, length(freqSamples));
                    for k = freqSamples
                        Phi_AZ = kron(squeeze(txGrid(:, :, k)).', W(:, :, k));
                        A_Theta = Phi_AZ*Psi_AZ;
                        run_MCs(k) = mutual_coherence(A_Theta);
                    end
                MCs(N_Theta_ind, U_ind, n) = mean(run_MCs, "all");
            end
        end
    end