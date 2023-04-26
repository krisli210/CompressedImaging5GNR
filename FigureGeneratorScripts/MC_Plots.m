close all

f = load('MCs.mat');

N_Theta_range = f.N_Theta_range;
U_range = f.U_range;

data = f.MCs;
data = data(:, :, 1:100); %oopsie 

figure; hold on;
for N_Theta_ind = 1:length(N_Theta_range)
    plot(U_range, mean(data(N_Theta_ind, :, :), 3), '-x');
end

legendString = string([repmat('$\delta_{\Theta} = ', length(N_Theta_range), 1), ...
    num2str(120./N_Theta_range.', 3), ...
    repmat('^{\circ}$ ', length(N_Theta_range), 1)]);
legend(legendString, 'Location', 'best', 'Interpreter','latex','Fontsize', 16);
xlabel('Users per Resource Block');
ylabel('$\mu(\mathbf{A}_{\Theta, k})$', 'Interpreter','latex')
title({'Mutual Coherence vs. $U$', ...
    ['$N_T = ', num2str(1) '; N, M = ', num2str(16), '$']}, 'Interpreter', 'latex');
grid on;