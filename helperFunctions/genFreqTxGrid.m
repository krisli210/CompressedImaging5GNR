function [txGrid, txAngleRecord] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook, Pt_W)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);
    
    sIndices = randi([0 MCS-1], [U, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); 
      
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this
    % txCodebook = txCodebook(:, 1:200);
    txAngleRecord = zeros(N_T, K, U); % Indices in H_RX
    for n_T = 1:N_T
        txAngles = randperm(size(txCodebook, 2), U); % Precoder changes per slot
        txAngleRecord(n_T, 1, :) = txAngles;
        F = 1./sqrt(U) .* txCodebook(:, txAngles); % M x U

        for nofdm = 1:Nofdm
            for k = 1:K
                if (~mod(k, 12+1))
                    txAngles = randperm(size(txCodebook, 2), U); % Precoder changes per 12 subcarriers
                    F = 1./sqrt(U) .* txCodebook(:, txAngles); % M x U
                end
                txAngleRecord(n_T, k, :) = txAngles;
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % U x 1
                unnormed_symbol_vector = F*s_slice;
                a = (Pt_W / M) * norm(unnormed_symbol_vector*unnormed_symbol_vector', 'fro');
                txGrid(:, startTimeIndex + nofdm, k) = a*unnormed_symbol_vector; % M x 1
            end
        end
    end
end