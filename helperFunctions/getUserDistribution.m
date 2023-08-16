function [theta_dist, theta_dist_cum, theta_dist_v, theta_dist_cum_v] = getUserDistribution(N_Theta, AzBins, distType)
    
    %% Generate custom bimodal
    if strcmp(distType, 'custom bimodal')
        n = 1e7;
        var = 40; % degrees
        mean1 = -45;
        mean2 = 45;
        bimodal_histogram = ( [sqrt(var) / 2 * randn(n/2, 1) - mean1
                                 sqrt(var) / 2 * randn(n/2, 1) - mean2 ] );
        
        % roundto = interp1(AzBins, AzBins, bimodal_histogram, 'nearest', 'extrap');
    
        % theta_dist = histogram(roundto, 'Normalization','pdf').Values;
        theta_dist = histc(bimodal_histogram, AzBins) / n;

    elseif (strcmp(distType, 'uniform'))
        theta_dist = ones(1, N_Theta) / N_Theta;
    end

    theta_dist_cum = cumsum(theta_dist);

    theta_dist_v = (max(theta_dist)-theta_dist) / sum(max(theta_dist)-theta_dist); 
    theta_dist_cum_v = cumsum(theta_dist_v);
end