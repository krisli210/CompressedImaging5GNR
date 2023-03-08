close all
clear 

rng(43);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 80; % number of resource blocks
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
    
    prm.freqSamples = 6; % number of subcarriers to sample in the rx
    % % % % % % % Grid/Target Construction
    prm.L = 2;
    prm.rMin = 20; prm.rMax = 100;

    limit_rangeRes = false;
    if (limit_rangeRes)
        prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*prm.K); % nominal range resolution
    else
        % Block for changing rangeRes to a non-nominal value
        prm.delta_R = 6;
    end
    prm.RangeBins = prm.rMin : prm.delta_R : prm.rMax;
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

% % % Transmit Signal Construction

    prm.NumUsers = 4;
    prm.N_T = 6; % number of time slots
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
        Y_tens(:, :, k) = rxGain .* (H_tens(:, :, k) * txGrid(:, :, k));
    end
    
% % % Receive Processing
    W = eye(prm.NumRxElements);
    % Az Cutting
    % Currently just grab any subcarrier - maybe has to do with not
    % dealing with scaled offset, dominated mostly by array structure
    kSubsample = randperm(prm.K, prm.freqSamples);
    x_AZ = txGrid(:, :, 1); % Space x Time - N x Nofdm * N_T
    y_AZ = Y_tens(:, :, 1);
    
    y_vec_AZ = reshape(y_AZ, [numel(y_AZ), 1]);

    Phi_AZ = kron(x_AZ.', W); % N*Nofdm*N_T x MN
    
    Psi_AZ = kr(H_TX, H_RX); % 
    Psi_AZ = Psi_AZ ./ norm(Psi_AZ);
    
    [z_hat_AZ, I_AZ] = solveCS_OMP(y_vec_AZ, Phi_AZ, Psi_AZ);

    % Construct ranging dic
    Psi_R = zeros(prm.K, prm.N_R);
    for r = 1:length(prm.RangeBins)
        tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
        Psi_R(:, r) = exp(-1j * 2*pi .* [0:prm.K-1]*prm.Delta_f * tau_r); % This is just the ifft of an identity
%         Psi_R(:, r) = Psi_R(:, r) ./ norm(Psi_R(:, r));
    end
    Psi_R = Psi_R ./ norm(Psi_R);
    Phi_R = eye(size(Psi_R, 1));

    % Range Cutting

    % Choose random antenna elements for now to test IFFT basis
    tx_rand = randi(prm.BsArraySize);
    rx_rand = randi(prm.RxArraySize);
    
    x_R = squeeze(txGrid(tx_rand, :, :)).'; % Freq x Time
    y_R = squeeze(Y_tens(rx_rand, :, :)).'; % 
    
%     weight_mags = ones(1, length(I_AZ)).';
%     y_R = steerTensor(Y_tens, H_RX, I_AZ, weight_mags).';
%     x_R = steerTensor(txGrid, H_TX, I_AZ, weight_mags).';
%     
    H_hat = mean(y_R ./ x_R, 2); % Average freq response over OFDM symbols
    [z_hat_R, I_R] = solveCS_OMP(H_hat, Phi_R, Psi_R);
%     z_hat_R = ifft(squeeze(mean(Y_tens(rx_rand, :, :) ./ txGrid(tx_rand, :, :), 2))); % Choose random antenna pairing, average over all OFDM symbols - this has higher interference compared to BF 
   
   
% % % END Receive Processing
    
figure; 
subplot(1, 2, 1); hold on; 
stem(-60:120/(prm.N_theta-1):60, abs(sum(RangeAzProfile, 1)));
stem(-60:120/(prm.N_theta-1):60, abs(z_hat_AZ), '--x');
xlabel('\theta');
ylabel('Magnitude');
title('Azimuth Profile')
legend({'True', 'Estimate'}, 'Location', 'best')  

subplot(1, 2, 2); hold on;
stem(prm.RangeBins, abs(sum(RangeAzProfile, 2)));
stem(prm.RangeBins, abs(z_hat_R), '--x');
xlabel('Range [m]');
ylabel('Magnitude');
title('Range Profile')
legend({'True', 'Estimate'}, 'Location', 'best') 

% Construct final image
L_hat = min(length(I_AZ), length(I_R));
[azMag_hat, azInd_hat] = maxk(z_hat_AZ, L_hat, 'ComparisonMethod','abs');
[rMag_hat, rInd_hat] = maxk(z_hat_R, L_hat, 'ComparisonMethod','abs');

RangeAzProfile_hat = zeros(size(RangeAzProfile));
for l = 1:L_hat
    RangeAzProfile_hat(rInd_hat(l), azInd_hat(l)) = azMag_hat(l);
end

figure; 
[h, c] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile_hat).^2).', ...
    'typerose', 'default', 'labelr', 'r [m]');
c.Label.String = 'Measured Reflection Power [dB]';