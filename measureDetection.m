function [detect, falseAlarm] = measureDetection(azProfile, z_hat)
    I = sort(find(azProfile));
    
    N = length(I); % num targets
    detect = zeros(1, N);
    for n = 1:N
        if (z_hat(I(n)) > .2) %arbitrary threshold 
            detect(n) = 1;
        end
    end

    z_hat(I) = [];
    falseAlarm = 0;
    if any(z_hat > .2) 
        falseAlarm = 1; 
    end
end