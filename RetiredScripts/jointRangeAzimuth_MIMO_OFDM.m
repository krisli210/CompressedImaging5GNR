close all
clear 

rng(42);


% % % carrier config

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120; % SCS in KHz
    
    prm.NRB = 30;
    prm.K = 12 * 30; %number of subcarriers = 12 subcarriers per RB * NRB - DOES NOT MATCH NFFT?
    carrier = nrCarrierConfig('SubcarrierSpacing', prm.Delta_f, 'NSizeGrid', prm.NRB);
    carrierInfo = nrOFDMInfo(carrier);
    prm.SampleRate = carrierInfo.SampleRate;
%     prm.K = carrierInfo.Nfft; % Number of subcarriers
    
% % % END carrier config

% % % % % Array and Channel Info

    prm.BsPos = [0; 0; 0];
    prm.BsArraySize = 16; %BS Dimension
    prm.NumBsElements = prod(prm.BsArraySize);
    prm.DeltaT = .5; % Element spacing normalized by wavelength
    prm.BsAZlim = [-60 60];
    prm.BsELlim = [-90 0];
    
    prm.RxPos = [0; 0; 0];
    prm.RxArraySize = 16;
    prm.DeltaR = prm.NumBsElements * .5; % Set for virtual array 
    prm.NumRxElements = prod(prm.RxArraySize);
    prm.RxAZlim = prm.BsAZlim;
    prm.RxELlim = [-90 0];
    
    %Arrays as uniform linear given in PA toolbox
    BsArray = phased.ULA(prm.NumBsElements, prm.DeltaT*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true), 'ArrayAxis', 'y');
    
    RxArray = phased.ULA(prm.NumRxElements, prm.DeltaR*prm.lam, 'Element', phased.IsotropicAntennaElement,'ArrayAxis', 'y');
    
    prm.N_theta = 64; %number of grid spots, i.e. dimension of the quantized azimuth profile
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;

    % Dictionary Construction
    H_TX = zeros(prm.NumBsElements, prm.N_theta);
    H_RX = zeros(prm.NumRxElements, prm.N_theta);
    for n = 1:prm.N_theta
        H_TX(:, n) = (1/sqrt(prm.NumBsElements)) * exp(1j * 2 * pi * prm.DeltaT * (0:prm.NumBsElements-1) * sind(prm.AzBins(n))).';
        H_RX(:, n) = (1/sqrt(prm.NumRxElements)) * exp(1j * 2 * pi * prm.DeltaR * (0:prm.NumRxElements-1) * sind(prm.AzBins(n))).';
%         H_TX(:, n) = collectPlaneWave(BsArray, 1, [prm.AzBins(n), 0].', prm.CenterFreq);
%         H_RX(:, n) = collectPlaneWave(RxArray, 1, [prm.AzBins(n), 0].', prm.CenterFreq);
    end
% % % % % END Array Info

% % % Transmit Signal Construction

    prm.NumUsers = 3;
    prm.N_T = 9; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 4; %modulation order

    txGrid = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX);

    % nrOFDMModulate takes in K x (Nofdm * N_T) x M outputs time-domain
    [txWaveform, info] = nrOFDMModulate(carrier, txGrid);
    
% % % END Transmit Signal Construction

% % % % % % % Target Construction
    prm.NumTargets = 1;
    prm.rMin = 40; prm.rMax = 100;

    [channel, azProfile, ScatPosCart, ScatPosPol] = genChannel(prm, BsArray, RxArray);
% % % % % % % END Target Construction

% % % Receive Processing

    [rxWaveform, chMat, tau] = channel(txWaveform);
    rxGrid = nrOFDMDemodulate(carrier, rxWaveform);
    
% Spatial Compression
    
%     for k = 1:prm.K
%         rxGrid(k, :, :) = exp(-1j * 2 * pi * tau * (prm.Delta_f * 1e3 * k)) .* rxGrid(k, :, :);
%     end
    
%     x = squeeze(mean(txGrid, 1)).';
%     freqAvg = squeeze(mean(rxGrid, 1));
    
    x = squeeze(txGrid(27, :, :)).';
    freqAvg = squeeze(rxGrid(27, :, :));
    W = eye(prm.NumRxElements);
    y_vec = reshape(freqAvg, [numel(freqAvg), 1]);

    Phi_az = kron(x.', W);
    Psi_az = kr(H_TX, H_RX);

    A = Phi_az * Psi_az;
    
    sensingDict = sensingDictionary('CustomDictionary', A);
    [z_hat, YI, I, R] = matchingPursuit(sensingDict, y_vec, maxIterations=10, Algorithm="OMP", maxerr={"L1", 1e-4});
    
    A_sub = Phi_az * Psi_az(:, I); % Solve magnitude posthence via direct linsolve against estimated support
    mags = linsolve(A_sub, y_vec);
    z_hat(I) = mags;
    
%     z_hat = linsolve(A, y_vec);
    figure; hold on; 
    stem(-60:120/(prm.N_theta-1):60, abs(azProfile)*max(abs(z_hat)), '--o');
    stem(-60:120/(prm.N_theta-1):60, abs(z_hat), '--x');
    xlabel('\theta');
    ylabel('Magnitude');
    title('Azimuth Profile')
    
    legend({'True', 'Estimate'})
    % END Spatial Compression

% % %

show_scene = false;
if (show_scene)
%     wT = BsSteeringVector(prm.CenterFreq,BsBeamAng(:,maxBlockRADAR));
%     wR = RxSteeringVector(prm.CenterFreq, RxBeamAng(:, maxBlockRADAR));
    txAngles = randperm(size(H_TX, 2), prm.NumUsers);
    F = 1./sqrt(prm.NumUsers) * H_TX(:, txAngles); % M x U
    wT = complex(randn(prm.NumBsElements, 1), randn(prm.NumBsElements, 1));
%     wR = complex(randn(64, 1), randn(64, 1));
    wR = .001*ones(prm.NumRxElements, 1);

    % Plot MIMO scenario with tx, rx, scatterers, and determined beams. Beam
    % patterns in this figure resemble the power patterns in linear scale.
    prmScene = struct();
    prmScene.TxArray = BsArray;
    prmScene.RxArray = RxArray;
    prmScene.TxArrayPos = prm.BsPos;
    prmScene.RxArrayPos = prm.RxPos;
    prmScene.ScatterersPos = ScatPosCart;
    prmScene.Lambda = prm.lam;
    prmScene.ArrayScaling = .1;     % To enlarge antenna arrays in the plot
    prmScene.MaxTxBeamLength = 20; % Maximum length of transmit beams in the plot
    prmScene.MaxRxBeamLength = 0; % Maximum length of receive beam in the plot
    hPlotSpatialMIMOScene(prmScene,F,wR);
%     if ~prm.ElevationSweep
%         view(2);
%     end
    view(2);
end
function [txGrid] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook)
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

function [channelRADAR, azProfile, ScatPosCart, ScatPosPol] = genChannel(prm, BsArray, RxArray)
    %Scatterer generation
    ScatCoeffs = complex(ones(1, prm.NumTargets), ones(1, prm.NumTargets)) ./ sqrt(2); %Unit reflector
    [ScatPosPol, ScatPosCart, azProfile] = generateDiscreteScatPositions(prm.NumTargets, prm.rMin, prm.rMax, prm.AzBins, ScatCoeffs);
    %End Scatterer generation
    
    % % % RADAR Channel Block
    channelRADAR = phased.ScatteringMIMOChannel;
    channelRADAR.PropagationSpeed = prm.PropagationSpeed;
    channelRADAR.CarrierFrequency = prm.CenterFreq;
    channelRADAR.SampleRate = prm.SampleRate;
    channelRADAR.SimulateDirectPath = false;
    channelRADAR.ChannelResponseOutputPort = true; %Gives delay as output
    channelRADAR.Polarization = 'None'; %Dual polarization requires defining scattering functions in both horiz and vert directions
    
    channelRADAR.TransmitArray = BsArray;
    channelRADAR.TransmitArrayPosition = prm.BsPos; 
    channelRADAR.ReceiveArray = RxArray;
    channelRADAR.ReceiveArrayPosition = prm.RxPos; 
    channelRADAR.ScattererSpecificationSource = 'Property';
    channelRADAR.ScattererPosition = ScatPosCart;
    channelRADAR.ScattererCoefficient = ScatCoeffs;
    channelRADAR.MaximumDelaySource = 'Auto';
   
end