close all
clear 

rng(42);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120; % SCS in KHz
    
    prm.K = 512; %number of subcarriers = 12 subcarriers per RB * NRB - DOES NOT MATCH NFFT?
    
% % % END OFDM Signal Params

% % % % % Array and Channel Info

    prm.BsPos = [0; 0; 0];
    prm.BsArraySize = 16; %BS Dimension
    prm.NumBsElements = prod(prm.BsArraySize);
    prm.DeltaTX = .5; % Element spacing normalized by wavelength
    prm.BsAZlim = [-60 60];
    prm.BsELlim = [-90 0];
    
    prm.RxPos = [0; 0; 0];
    prm.RxArraySize = 15;
    prm.DeltaRX = prm.NumBsElements * .5; % Set for virtual array 
    prm.NumRxElements = prod(prm.RxArraySize);
    prm.RxAZlim = prm.BsAZlim;
    prm.RxELlim = [-90 0];
    
    %Arrays as uniform linear given in PA toolbox
    BsArray = phased.ULA(prm.NumBsElements, prm.DeltaTX*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true), 'ArrayAxis', 'y');
    
    RxArray = phased.ULA(prm.NumRxElements, prm.DeltaRX*prm.lam, 'Element', phased.IsotropicAntennaElement,'ArrayAxis', 'y');
    
    prm.N_theta = 64; %number of grid spots, i.e. dimension of the quantized azimuth profile
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    
    M = 1; %ifft oversampling factor
    prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*1e3*prm.K); % nominal range resolution
    prm.N_R = 64;
    prm.RangeBins = prm.delta_R:prm.delta_R:(prm.N_R-1)*prm.delta_R;
    maxRange = prm.RangeBins(end);

    % Dictionary Construction
    H_TX = zeros(prm.NumBsElements, prm.N_theta);
    H_RX = zeros(prm.NumRxElements, prm.N_theta);
    for n = 1:prm.N_theta
%         H_TX(:, n) = (1/sqrt(prm.NumBsElements)) * exp(1j * 2 * pi * prm.DeltaT * (0:prm.NumBsElements-1) * sind(prm.AzBins(n))).';
%         H_RX(:, n) = (1/sqrt(prm.NumRxElements)) * exp(1j * 2 * pi * prm.DeltaR * (0:prm.NumRxElements-1) * sind(prm.AzBins(n))).';
        H_TX(:, n) = collectPlaneWave(BsArray, 1, [prm.AzBins(n), 0].', prm.CenterFreq);
        H_RX(:, n) = collectPlaneWave(RxArray, 1, [prm.AzBins(n), 0].', prm.CenterFreq);
    end
% % % % % END Array Info

% % % Transmit Signal Construction

    prm.NumUsers = 3;
    prm.N_T = 9; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 16; %modulation order
    
    %Generates baseband frequency-domain signal across freq-time-space
    txGrid = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % K x (Nofdm * N_T) x M
    
% % % END Transmit Signal Construction

% % % % % % % Target Construction
    prm.L = 1;
    prm.rMin = 40; prm.rMax = 100;

    [H_tens, RangeAzProfile, ScatPosPol] = genGridChannel(prm);
% % % % % % % END Target Construction

% % % Receive Processing
    
% % %

function [txGrid] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is K x (Nofdm * N_T) x M
    txGrid = zeros(K, Nofdm * N_T, M);

    sIndices = randi([0 MCS-1], [K, Nofdm * N_T, U]);     % per user symbols given as K x Nofdm * N_T x U
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); 
    
    % precoded transmit symbols given as K x Nofdm * N_T x M

    % Loop over subcarriers and slots because idk how to tensorize this
    for k = 1:K % FIXME: actually beamforming per RB 
        for n_T = 1:N_T
            txAngles = randperm(size(txCodebook, 2), U);
            F = 1./sqrt(U) * txCodebook(:, txAngles); % M x U
            for nofdm = 1:Nofdm
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(k, startTimeIndex + nofdm, :)); % U x 1
                txGrid(k, startTimeIndex + nofdm, :) = F*s_slice;
            end
        end
    end
end

function [H_tens, RangeAzProfile, ScatPosPol] = genGridChannel(prm)
    % H _tens output as N x M x K
    % RangeAzProfile as N_R x N_theta
    % ScatPosPol as 3 x L
    azInd = randperm(prm.N_theta, prm.L);
    rangeInd = randperm(prm.N_R, prm.L);
    
    azValues = prm.AzBins(azInd);
    rangeValues = prm.RangeBins(rangeInd);
    ScatPosPol = [rangeValues; azValues; zeros(1, prm.L)];
    ScatCoeff = ones(1, 3) .* complex(1, 1) ./ sqrt(2); %Unit reflectors

    RangeAzProfile = zeros(prm.N_R, prm.N_theta);
    H_tens = zeros(prm.NumRxElements, prm.NumBsElements, prm.K);
    for l = 1:prm.L
        RangeAzProfile(rangeInd(l), azInd(l)) = ScatCoeff(l);
        tau_r = 2*rangeValues(l) / prm.PropagationSpeed;
        tau_m = prm.DeltaTX * (0:prm.NumBsElements-1)*sind(azValues(l)); %Check this
        tau_n = prm.DeltaRX * (0:prm.NumRxElements-1)*sind(azValues(l));
        tau_n_m = zeros(prm.NumRxElements, prm.NumBsElements);

        for n = 1:prm.NumRxElements
            for m = 1:prm.NumBsElements
                tau_n_m(n,m) = tau_n(n) + tau_m(m);
            end
        end

        tau_total = tau_r + tau_n_m;
        PL = (4*pi*rangeValues(l)/prm.lam)^-2;
        for k = 1:prm.K
            H_tens(:, :, k) = H_tens(:, :, k) + PL*ScatCoeff(l)*exp(-1j * 2*pi * (k*prm.Delta_f*1e3) .* tau_total);
        end
    end
end