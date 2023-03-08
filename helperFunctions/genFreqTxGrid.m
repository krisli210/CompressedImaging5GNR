function [txGrid, txGridLB] = genFreqTxGrid(M, U, MCS, N_T, Nofdm, K, txCodebook, varargin)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);
    
    if (nargin == 8)
        azInd = varargin{1};
        L = length(azInd);
        
        txGridLB = zeros(size(txGrid)); % We use as a transmit signal sent only to the scattering centers - use as a lower bound on error
    end

    sIndices = randi([0 MCS-1], [U, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); 
    
    if (nargin == 8)
        sIndicesLB = randi([0 MCS-1], [L, Nofdm * N_T, K]);
        sLB = qammod(sIndicesLB, MCS, 'UnitAveragePower', true); 
        FLB = 1./sqrt(L) * txCodebook(:, azInd); % M x L - only set one as it is optim
    end
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this

    for n_T = 1:N_T
        txAngles = randperm(size(txCodebook, 2), U); % Precoder changes per slot 
        F = 1./sqrt(U) * txCodebook(:, txAngles); % M x U
        for nofdm = 1:Nofdm
            for k = 1:K
                if (~mod(k, 12+1))
                    txAngles = randperm(size(txCodebook, 2), U); % Precoder changes per 12 subcarriers
                    F = 1./sqrt(U) * txCodebook(:, txAngles); % M x U
                end
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % U x 1
                txGrid(:, startTimeIndex + nofdm, k) = F*s_slice;
                
                if (nargin == 8)
                    s_sliceLB = squeeze(sLB(:, startTimeIndex + nofdm, k)); % L x 1
                    txGridLB(:, startTimeIndex + nofdm, k) = FLB*s_sliceLB;
                end
            end
        end
    end
end