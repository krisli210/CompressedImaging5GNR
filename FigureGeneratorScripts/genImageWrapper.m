rng(42);

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

% % % % % % % Target Construction
%     prm.N_theta = prm.BsArraySize * prm.RxArraySize; %number of grid spots, i.e. dimension of the quantized azimuth profile
prm.N_theta = 256;
thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;

% % % % % % % Grid/Target Construction
prm.L = 50;
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

nIters = 1000;
SNR_range = -20:5:20;
container = zeros(3, length(SNR_range), nIters);
for SNR_ind = 1:length(SNR_range)
    prm.SNR_dB = SNR_range(SNR_ind);
      
    for i = 1:nIters
        container(1, SNR_ind, i) = genImageWrapped_SS_randW(prm, H_TX, Psi_AZ, Psi_R, 'eye');
        container(2, SNR_ind, i) = genImageWrapped_SS_randW(prm, H_TX, Psi_AZ, Psi_R, 'ffteye');
        container(3, SNR_ind, i) = genImageWrapped_SS_randW(prm, H_TX, Psi_AZ, Psi_R, 'rand');
    end
end
