function [txGrid, txGridLB] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook, azInd)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);
    txGridLB = zeros(size(txGrid)); % We use as a transmit signal sent only to the scattering centers - use as a lower bound on error
    L = length(azInd);

    sIndices = randi([0 MCS-1], [U, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); 
    
    sIndicesLB = randi([0 MCS-1], [L, Nofdm * N_T, K]);
    sLB = qammod(sIndicesLB, MCS, 'UnitAveragePower', true); 
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this
    for k = 1:K % FIXME: actually beamforming per RB 
        for n_T = 1:N_T
            txAngles = randperm(size(txCodebook, 2), U);
            F = 1./sqrt(U) * txCodebook(:, txAngles); % M x U - Change per frame (14 OFDM symbols) and per subcarrier (although should be per 12 subcarriers)
            FLB = 1./sqrt(L) * txCodebook(:, azInd); % M x L
            for nofdm = 1:Nofdm
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % U x 1
                s_sliceLB = squeeze(sLB(:, startTimeIndex + nofdm, k)); % L x 1

                txGrid(:, startTimeIndex + nofdm, k) = F*s_slice;
                txGridLB(:, startTimeIndex + nofdm, k) = FLB*s_sliceLB;
            end
        end
    end
end