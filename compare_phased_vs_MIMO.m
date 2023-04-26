close all
f = load("compare_phased_vs_MIMO_U_SNR_N_T_1.mat");

data = squeeze(f.wholeData);
delta_Theta = 120/f.prm.N_theta;

U_range = f.U_range;
SNR_range = f.SNR_dB_range;

figure; hold on;
for U_ind = 1:length(U_range)
    plot(SNR_range(1:9), 10*log10(mean(data(U_ind, 1:9, :), 3)), '-x');
end
grid on;

legendString = ([repmat('$U = ', [3 1]), num2str(U_range.'), repmat('$', [3 1])]);
legend(legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('Receive SNR [dB]')
ylabel('Mean NSE [dB]');
title({'Normalized Squared Error vs. $SNR, U$', ...
    ['$L = ', num2str(50) '; \delta_\Theta = ', num2str(delta_Theta, 3), '^{\circ}$']}, 'Interpreter', 'latex');
