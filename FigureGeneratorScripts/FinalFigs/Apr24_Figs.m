close

f = load('NSEs_L_U_noisy.mat');
U_range = f.U_range;
L_range = f.L_range;
delta_Theta = 120/f.N_Theta_range;

figure; hold on;
data = squeeze(f.wholeData);
for U_ind = 1:length(U_range)
    plot(L_range, 10*log10(mean(data(:, U_ind, :), 3)), '-x');
end
legendString = ([repmat('$U = ', [3 1]), num2str([1 3 5].'), repmat(' $', [3 1])]);
legend(legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
xlabel('$L$', 'Interpreter','latex')
ylabel('Mean NSE [dB]');
title({'Normalized Squared Error vs. $L, U$', ...
    ['$N_T = ', num2str(1) '; \delta_\Theta = ', num2str(delta_Theta, 3), '^{\circ}$']}, 'Interpreter', 'latex');