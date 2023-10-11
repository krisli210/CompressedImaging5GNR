close all
clear

% MNSE vs. N_T 
load("NSEs_SNRs_SS_Comparison_vary_N_theta.mat");
NSEs_SNRs = 10*log10(NSEs_SNRs(1:length(SNR_range), :, :));
NSE_SNRs_SS = 10*log10(NSEs_SNRs_SS(1:length(SNR_range), :, :));

N_theta_range = N_theta_range(1:end-1);
col = ['r', 'g', 'b', 'c'];
figure; hold on; 
for i = 1:length(N_theta_range)
    plot(SNR_range, mean(NSEs_SNRs(:, i, :), 3), '-o', 'color', col(i));
    plot(SNR_range, mean(NSE_SNRs_SS(:, i, :), 3), '--x', 'color', col(i));
end

for i = 1:length(N_theta_range)
    qw{i} = plot(nan, 'color', col(i));
end
qw{i+1} = plot(SNR_range, 100*ones(1,length(SNR_range)), '-o', 'color', [0 0 0]);
qw{i+2} = plot(SNR_range, 100*ones(1,length(SNR_range)), '--x', 'color', [0 0 0]);
grid on;
xlabel('Receive SNR [dB]');
ylabel('Mean NSE [dB]');
ylim([-40 5]);
title({'Effect of Frequency Subsampling', '$N_T = 1; U=3$'}, 'Interpreter', 'latex');
legendString = string([repmat('$\delta_{\Theta} =  ', length(N_theta_range), 1), ...
    num2str(120./N_theta_range.', 3), repmat('^\circ$ ', length(N_theta_range), 1)]);
legendString = [legendString.', string('$K^\prime = K$'), string('$K^\prime = N_\mathcal{R}$')].';
legend([qw{:}], legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 12);