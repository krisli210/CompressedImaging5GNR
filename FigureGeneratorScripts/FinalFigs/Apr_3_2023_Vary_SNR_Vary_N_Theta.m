close all
clear

% MNSE vs. N_T 
load("NSEs_SNRs_SS_Comparison_vary_N_theta.mat");
NSEs_SNRs = NSEs_SNRs(1:length(SNR_range), :, :);
NSE_SNRs_SS = NSEs_SNRs_SS(1:length(SNR_range), :, :);

col = ['r', 'g', 'b', 'c'];
figure; hold on; 
for i = 1:length(N_theta_range)
    plot(SNR_range, mean(NSEs_SNRs(:, i, :), 3), '-o', 'color', col(i));
    plot(SNR_range, mean(NSE_SNRs_SS(:, i, :), 3), '--x', 'color', col(i));
end

for i = 1:length(N_theta_range)
    qw{i} = plot(nan, 'color', col(i));
end
grid on;
xlabel('Receive SNR [dB]');
ylabel('Normalized Squared Error');
title('Effect of Receive SNR on Image Accuracy')
legendString = string([repmat('$\delta_{\Theta} =  ', length(N_theta_range), 1), ...
    num2str(120./N_theta_range.', 2), repmat('^\circ$ ', length(N_theta_range), 1)]);
legend([qw{:}], legendString, 'location', 'best', 'Interpreter', 'latex');