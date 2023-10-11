rng(42);
warning('off')
clear

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

% prm.Pt = 20; %dBm 
% prm.Pr = 20; %dBm
% 
% prm.No = -174; %dBm/Hz
prm.SNR_dB = -10; % Defined as the receive signal power / noise variance
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

thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
% Fixed Parameters
prm.N_T = 1;
prm.NumUsers = 3;
prm.SNR_dB = 10;

prm.rMin = 20; prm.rMax = 70;
max_FSPL_dB = 10*log10((4*pi*prm.rMax/prm.lam)^-2);
%

% Ranged Parameters
% L_range = [10:10:100];
L_range = 50;

N_Theta_range = [256];
% N_Theta_range = [64 128 256];

delta_R_range = 2;
N_R_range = zeros(1, length(delta_R_range));

N_T_range = [1:8];
% N_T_range = 1;

U_range = [1:2:5];
% U_range = 3;
nIters = 100;

NSEs_NT_NTheta_SS = zeros(length(N_Theta_range), length(delta_R_range), length(L_range), length(N_T_range), length(U_range));
wholeData = zeros([size(NSEs_NT_NTheta_SS), nIters]);

for N_Theta_ind = 1:length(N_Theta_range)

    prm.N_theta = N_Theta_range(N_Theta_ind); % % % % 

    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    H_TX = zeros(prm.NumBsElements, prm.N_theta);
    H_RX = zeros(prm.NumRxElements, prm.N_theta);
    for n = 1:prm.N_theta
        H_TX(:, n) = exp(-1j * 2 * pi * prm.DeltaTX * (0:prm.NumBsElements-1) * sind(prm.AzBins(n))).';
        H_RX(:, n) = exp(-1j * 2 * pi * prm.DeltaRX * (0:prm.NumRxElements-1) * sind(prm.AzBins(n))).';
    end

    Psi_AZ = kr(H_TX, H_RX); %

    for delta_R_ind = 1:length(delta_R_range)

        prm.delta_R = delta_R_range(delta_R_ind); % % % 

        prm.WholeRange = 0:prm.delta_R:(prm.K-1)*prm.delta_R;
        minIndex = find(prm.WholeRange < prm.rMin, 1, 'last')+1;
        maxIndex = find(prm.WholeRange > prm.rMax, 1, 'first')-1;
    
        prm.RangeBins = prm.WholeRange(minIndex:maxIndex);
        prm.N_R = numel(prm.RangeBins);
        N_R_range(delta_R_ind) = prm.N_R;

        Psi_R = zeros(prm.K, prm.N_R);
        for r = 1:length(prm.RangeBins)
            tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
            Psi_R(:, r) = exp(-1j * 2*pi .* [1:prm.K]*prm.Delta_f * tau_r); % This is just the ifft of an identity
        end

        for L_ind = 1:length(L_range)

            prm.L = L_range(L_ind); % % % 
            for N_T_ind = 1:length(N_T_range)
                prm.N_T = N_T_range(N_T_ind);
                for U_ind = 1:length(U_range)
                    prm.NumUsers = U_range(U_ind);
                    for i = 1:nIters
                        wholeData(N_Theta_ind, delta_R_ind, L_ind, N_T_ind, U_ind, i) = genImageWrapped_SS(prm, H_TX, Psi_AZ, Psi_R);
                    end
                    [NSEs_NT_NTheta_SS(N_Theta_ind, delta_R_ind, L_ind, N_T_ind, U_ind)] = mean(wholeData(N_Theta_ind, delta_R_ind, L_ind, N_T_ind, U_ind, :), 'all');
                end
            end
        end
    end
end

save('NSEs_NT_NTheta_noiseless', "NSEs_NT_NTheta_SS", 'wholeData', 'N_Theta_range', 'delta_R_range', 'L_range', 'N_T_range', 'U_range', 'prm');