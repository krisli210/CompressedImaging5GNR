function [txGrid, powerLog, userAngleIndLog] = genAngleDefFreqTxGrid(M, U, V, MCS, N_T, Nofdm, K, txCodebook, Pt_W, theta_dist_cum)
    N_Theta = size(txCodebook, 2);
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);
    powerLog = [];
    userAngleIndLog = [];

    sIndices = randi([0 MCS-1], [N_Theta, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); % this should take care of E[ss^H] = I_U / U ?
      
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this
    for n_T = 1:N_T
        
        [p, userAngleInds] = getUserPowerVector(U, theta_dist_cum, N_Theta, Pt_W);
        powerLog = [powerLog p(userAngleInds)];
        userAngleIndLog = [userAngleIndLog userAngleInds];
        % F = 1./sqrt(U) .* txCodebook(:, userAngles); % M x U

        for nofdm = 1:Nofdm
            for k = 1:K
                if (~mod(k, 12+1))
                    [p, userAngleInds] = getUserPowerVector(U, theta_dist_cum, N_Theta, Pt_W);
                    powerLog = [powerLog p(userAngleInds)];
                    userAngleIndLog = [userAngleIndLog userAngleInds];
                end
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % N_Theta x 1

                % s_slice_expand = zeros(N_Theta, 1);
                % s_slice_expand(userAngleInds) = s_slice;

                txGrid(:, startTimeIndex + nofdm, k) = txCodebook * sqrt(diag(p)) * s_slice;
            end
        end
    end
end

function [p, userAngleInds] = getUserPowerVector(U, theta_dist_cum, N_Theta, Pt_W)
    vals = rand(U, 1);
    randfun = @(r) find(r <= theta_dist_cum, 1, 'first');
    userAngleInds = arrayfun(randfun, vals);
    
    % power matrix definition according to non-zero angle bins by users
    p = zeros(N_Theta, 1);
    p(userAngleInds) = Pt_W / numel(userAngleInds); % even power distribution among real users
    
    % p = Pt_W / N_Theta * ones(N_Theta, 1);
end