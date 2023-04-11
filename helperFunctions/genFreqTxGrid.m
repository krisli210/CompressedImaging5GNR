function [txGrid] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);
    
    sIndices = randi([0 MCS-1], [U, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); 
      
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this

    for n_T = 1:N_T
        txAngles = randperm(size(txCodebook, 2), U); % Precoder changes per slot
        F = 1./sqrt(U) .* txCodebook(:, txAngles); % M x U

        for nofdm = 1:Nofdm
            for k = 1:K
                if (~mod(k, 12+1))
                    txAngles = randperm(size(txCodebook, 2), U); % Precoder changes per 12 subcarriers
                    F = 1./sqrt(U) .* txCodebook(:, txAngles); % M x U
                end
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % U x 1
                txGrid(:, startTimeIndex + nofdm, k) = F*s_slice;
            end
        end
    end
end