close all
clear 

rng(44);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 40; % number of resource blocks
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
    prm.L = 4;
    prm.rMin = 20; prm.rMax = 100;

    limit_rangeRes = true;
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
    
    Psi_AZ = kr(H_TX, H_RX); % 
    Psi_AZ = Psi_AZ ./ norm(Psi_AZ);

    Psi_R = zeros(prm.K, prm.N_R);
    for r = 1:length(prm.RangeBins)
        tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
        Psi_R(:, r) = exp(-1j * 2*pi .* [0:prm.K-1]*prm.Delta_f*tau_r); % This is just the ifft of an identity
%         Psi_R(:, r) = Psi_R(:, r) ./ norm(Psi_R(:, r));
    end
    Psi_R = Psi_R ./ norm(Psi_R);
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
    RangeAzProfile_hat = sequentialOMP(txGrid, Y_tens, Psi_AZ, Psi_R, H_TX, H_RX, prm);
    
% % % END Receive Processing
    
 

% Construct final image

% RangeAzProfile_hat = zeros(size(RangeAzProfile));
% for l = 1:L_hat
%     RangeAzProfile_hat(rInd_hat(l), azInd_hat(l)) = azMag_hat(l);
% end

figure;
subplot(1, 2, 1);
[h, c] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile).^2).', ...
    'typerose', 'default', 'labelr', 'r [m]');
c.Label.String = 'Measured Reflection Power [dB]';

subplot(1, 2, 2);
[h, c] = polarPcolor(prm.RangeBins, prm.AzBins, 10*log10(abs(RangeAzProfile_hat).^2).', ...
    'typerose', 'default', 'labelr', 'r [m]');
c.Label.String = 'Measured Reflection Power [dB]';

sgtitle('True (L) vs. Estimate (R)'); 