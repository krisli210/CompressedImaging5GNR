function [theta_dist, theta_dist_cum] = getUserDistribution(N_Theta, AzBins, distType)
    
    %% Generate custom bimodal
    if strcmp(distType, 'custom bimodal')
        n = 1e7;
        var = 40; % degrees
        mean1 = -30;
        mean2 = 30;
        bimodal_histogram = ( [sqrt(var) / 2 * randn(n/2, 1) - mean1
                                 sqrt(var) / 2 * randn(n/2, 1) - mean2 ] );
        
        % roundto = interp1(AzBins, AzBins, bimodal_histogram, 'nearest', 'extrap');
    
        % theta_dist = histogram(roundto, 'Normalization','pdf').Values;
        theta_dist = histc(bimodal_histogram, AzBins) / n;

    elseif (strcmp(distType, 'uniform'))
        theta_dist = ones(1, N_Theta) / N_Theta;
    end

    theta_dist_cum = cumsum(theta_dist);
end