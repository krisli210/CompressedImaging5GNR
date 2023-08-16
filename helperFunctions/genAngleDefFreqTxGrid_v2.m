function [txGrid, powerLog, userAngleIndLog] = genAngleDefFreqTxGrid_v2(M, U, V, alpha, MCS, N_T, Nofdm, K, txCodebook, Pt_W, theta_dist_cum, theta_dist_cum_v)
    N_Theta = size(txCodebook, 2);
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(M, Nofdm * N_T, K);
    powerLog = [];
    userAngleIndLog = [];

    sIndices = randi([0 MCS-1], [U+V, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); % this should take care of E[ss^H] = I_U / U ?
      
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this
    for n_T = 1:N_T
        
        [p, userAngleInds, userAngleInds_v] = getUserPowerVector(U, V, alpha, theta_dist_cum, theta_dist_cum_v, Pt_W);
        powerLog = [powerLog p(1:U)];
        userAngleIndLog = [userAngleIndLog userAngleInds];
        F = txCodebook(:, [userAngleInds.' userAngleInds_v.']);    

        for nofdm = 1:Nofdm
            for k = 1:K
                if (~mod(k, 12+1))
                    [p, userAngleInds, userAngleInds_v] = getUserPowerVector(U, V, alpha, theta_dist_cum, theta_dist_cum_v, Pt_W);
                    powerLog = [powerLog p(1:U)];
                    userAngleIndLog = [userAngleIndLog userAngleInds];
                end
                startTimeIndex = (Nofdm * (n_T - 1));
                s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % N_Theta x 1

                F = txCodebook(:, [userAngleInds.' userAngleInds_v.']);
                txGrid(:, startTimeIndex + nofdm, k) = F * sqrt(diag(p)) * s_slice;
            end
        end
    end
end

function [p, userAngleInds, userAngleInds_v] = getUserPowerVector(U, V, alpha, theta_dist_cum, theta_dist_cum_v, Pt_W)
    vals = rand(U, 1); % 
    userAngleInds = zeros(U, 1);
    for u = 1:U
        userAngleInds(u) = find(vals(u) <= theta_dist_cum, 1, 'first');
    end
    
    vals_virtual = rand(V, 1);
    userAngleInds_v = zeros(V, 1);
    for v = 1:V
        userAngleInds_v(v) = find(vals_virtual(v) <= theta_dist_cum_v, 1, 'first');
    end
    % power matrix definition according to non-zero angle bins by users
    p = zeros(U+V, 1);
    p(1:U) = alpha * Pt_W / numel(userAngleInds); % even power distribution among real users
    p(U+1:U+V) = (1-alpha) * Pt_W / numel(userAngleInds_v);

    % p = Pt_W / N_Theta * ones(N_Theta, 1);
end