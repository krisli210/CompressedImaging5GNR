close all
clear 

rng(42);


% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.K = 512;
    
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
    prm.N_theta = 128;
    thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;
    
    % % % % % % % Grid/Target Construction
    prm.L = 3;
    prm.rMin = 20; prm.rMax = 150;

    limit_rangeRes = true;
    if (limit_rangeRes)
        prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*prm.K); % nominal range resolution
    else
        % Block for changing rangeRes to a non-nominal value
        prm.delta_R = 1;
        
    end
    prm.RangeBins = prm.rMin : prm.delta_R : prm.rMax;
    prm.N_R = numel(prm.RangeBins);

    [H_tens, RangeAzProfile, ScatPosPol] = genGridChannel(prm);
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
    prm.N_T = 30; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 16; %modulation order
    
    %Generates baseband frequency-domain signal across freq-time-space
    txGrid = genFreqTxGrid(prm.NumBsElements, prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K, H_TX); % (Nofdm * N_T) x M x K
    
% % % END Transmit Signal Construction

    Y_tens = zeros(prm.NumRxElements, prm.Nofdm * prm.N_T, prm.K);
    for k = 1:prm.K
        Y_tens(:, :, k) = (H_tens(:, :, k) * txGrid(:, :, k));
    end
    
% % % Receive Processing
    W = eye(prm.NumRxElements);
    % Az Cutting

    % Currently just grab any subcarrier - maybe has to do with not
    % dealing with scaled offset, dominated mostly by array structure
    x_AZ = squeeze(txGrid(:, :, 1)); % Space x Time - N x Nofdm * N_T
    y_AZ = squeeze(Y_tens(:, :, 1));
    
%     x_AZ = squeeze( mean ( txGrid(:, :, 1), 2)); % Broken for some reason
%     y_AZ = squeeze( mean ( Y_tens(:, :, 1), 2));
    
    y_vec_AZ = reshape(y_AZ, [numel(y_AZ), 1]);

    Phi_AZ = kron(x_AZ.', W); % N*Nofdm*N_T x MN
    Psi_AZ = kr(H_TX, H_RX); % 
    Psi_AZ = Psi_AZ ./ norm(Psi_AZ);
    
    % Construct ranging dic
    Psi_R = zeros(prm.K, prm.N_R);
    for r = 1:length(prm.RangeBins)
        tau_r = 2 * prm.RangeBins(r) / prm.PropagationSpeed;
        Psi_R(:, r) = exp(-1j * 2*pi .* [0:prm.K-1]*prm.Delta_f * tau_r); % This is just the ifft of an identity
        Psi_R(:, r) = Psi_R(:, r) ./ norm(Psi_R(:, r));
    end
    Phi_R = eye(size(Psi_R, 1));

    % Range Cutting
    
% % %

    % Choose random antenna elements for now to test IFFT basis
    tx_rand = randi(prm.BsArraySize);
    rx_rand = randi(prm.RxArraySize);
    
    x_R = squeeze(txGrid(tx_rand, :, :)).'; % Freq x Time
    y_R = squeeze(Y_tens(rx_rand, :, :)).'; % 
    
    y_R_vec = reshape(y_R.', [numel(y_R), 1]);
%     y_R = steerTensor(Y_tens, H_RX, I_AZ, z_hat_AZ(I_AZ)).';
%     x_R = steerTensor(txGrid, H_TX, I_AZ, z_hat_AZ(I_AZ)).';
    
    H_hat = mean(y_R ./ x_R, 2); % Average freq response over OFDM symbols
    
    [z_hat_AZ, I_AZ] = solveCS_OMP(y_vec_AZ, Phi_AZ, Psi_AZ);
    [z_hat_R, I_R] = solveCS_OMP(H_hat, Phi_R, Psi_R);
%     z_hat_R = linsolve(Psi_R, H_hat);
%     z_hat_R = ifft(squeeze(mean(Y_tens(rx_rand, :, :) ./ txGrid(tx_rand, :, :), 2))); % Choose random antenna pairing, average over all OFDM symbols - this has higher interference compared to BF 
    
    figure; 
    subplot(1, 2, 1); hold on; 
    stem(-60:120/(prm.N_theta-1):60, abs(sum(RangeAzProfile, 1)));
    stem(-60:120/(prm.N_theta-1):60, abs(z_hat_AZ), '--x');
    xlabel('\theta');
    ylabel('Magnitude');
    title('Azimuth Profile')
    legend({'True', 'Estimate'}, 'Location', 'south')  
    
    subplot(1, 2, 2); hold on;
    stem(prm.RangeBins, abs(sum(RangeAzProfile, 2)));
    stem(prm.RangeBins, abs(z_hat_R), '--x');
    xlabel('Range [m]');
    ylabel('Magnitude');
    title('Range Profile')
    legend({'True', 'Estimate'}, 'Location', 'south') 
% % % 
function [txGrid] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);

    sIndices = randi([0 MCS-1], [U, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); 
    
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this
    for k = 1:K % FIXME: actually beamforming per RB 
        for n_T = 1:N_T
            txAngles = randperm(size(txCodebook, 2), U);
            F = 1./sqrt(U) * txCodebook(:, txAngles); % M x U - Change per frame (14 OFDM symbols) and per subcarrier (although should be per 12 subcarriers)
            for nofdm = 1:Nofdm
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % U x 1
                txGrid(:, startTimeIndex + nofdm, k) = F*s_slice;
            end
        end
    end
end

function [H_tens, RangeAzProfile, ScatPosPol] = genGridChannel(prm)
    % H _tens output as N x M x K
    % RangeAzProfile as N_R x N_theta
    % ScatPosPol as 3 x L
    minRangeBin = find(prm.RangeBins <= prm.rMin, 1, 'last') + 1;
    maxRangeBin = find(prm.RangeBins >= prm.rMax, 1, 'first');
%     rangeInd = randperm(maxRangeBin, prm.L);
    rangeInd = randi([minRangeBin, maxRangeBin], [1, prm.L]);

    azInd = randperm(prm.N_theta, prm.L);
    rangeInd = randperm(prm.N_R, prm.L);

    azValues = prm.AzBins(azInd);
%     azValues = zeros(1, prm.L); % Fix az to 0 to test ranging
    rangeValues = prm.RangeBins(rangeInd);
%     rangeValues = ones(1, prm.L); % Fix range to 0 to test spatial compression
    ScatPosPol = [rangeValues; azValues; zeros(1, prm.L)];
    ScatCoeff = ones(1, prm.L) .* complex(1, 1) ./ sqrt(2); %Unit reflectors

    RangeAzProfile = zeros(prm.N_R, prm.N_theta);
    H_tens = zeros(prm.NumRxElements, prm.NumBsElements, prm.K);
    for k = 0:prm.K-1
        for l = 1:prm.L
                tau_r = 2*rangeValues(l) / prm.PropagationSpeed;
                tau_m = prm.DeltaTX * (0:prm.NumBsElements-1)*sind(azValues(l)); %Check this
                tau_n = prm.DeltaRX * (0:prm.NumRxElements-1)*sind(azValues(l));
                tau_n_m = zeros(prm.NumRxElements, prm.NumBsElements);
        
                for n = 1:prm.NumRxElements
                    for m = 1:prm.NumBsElements
                        tau_n_m(n,m) = tau_n(n) + tau_m(m);
                    end
                end
        
                PL = (4*pi*rangeValues(l)/prm.lam)^-2;
                RangeAzProfile(rangeInd(l), azInd(l)) = PL*ScatCoeff(l);
                H_tens(:, :, k+1) = H_tens(:, :, k+1) + PL*ScatCoeff(l)*exp(-1j * 2*pi * ( (k*prm.Delta_f*tau_r) + tau_n_m) ); % Is the az-induced delay scaled by freq? - NO
        end
    end
    H_tens = H_tens ./ prm.NumBsElements; %Power norm
end