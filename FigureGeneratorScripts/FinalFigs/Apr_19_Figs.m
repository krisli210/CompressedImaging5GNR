close all

f = load("NSEs_L_HighSNR_SS.mat");

L_range = f.L_range;
data = squeeze(f.NSEs_L_Res_SS);
delta_Theta  = 120 / f.N_Theta_range;
delta_R = f.delta_R_range;

figure;
plot(L_range, 10*log10(data), '-x');
xlabel('L');
ylabel('Mean NSE [dB]')
title({'Normalized Squared Error vs. L', ...
    ['$\delta_\Theta = ', num2str(delta_Theta, 3), ...
    '^{\circ} ; \delta_\mathcal{R} = ', num2str(delta_R), '[m]$']}, 'Interpreter', 'latex');
grid on;

f = load("NSEs_N_Theta_SS.mat");
N_Theta_range = f.N_Theta_range;

data = squeeze(f.NSEs_L_Res_SS);
L = f.prm.L;
figure;
plot(N_Theta_range, 10*log10(data), '-o');
xlabel({'$N_\Theta$'}, 'Interpreter', 'latex');
ylabel('Mean NSE [dB]');
title({'Normalized Squared Error vs. $N_\Theta$', ['$L = ', num2str(L) '$']}, 'Interpreter', 'latex');
grid on;

f = load('NSEs_N_T_U_noisy.mat');
U_range = f.U_range;
N_T_range = f.N_T_range;

figure; hold on;
data = squeeze(f.NSEs_L_U_SS);
for U_ind = 1:length(U_range)
    plot(N_T_range, 10*log10(data(:, U_ind)), '-x');
end
legendString = ([repmat('$U = ', [3 1]), num2str([1 3 5].'), repmat(' $', [3 1])]);
legend(legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
xlabel('$N_T$', 'Interpreter','latex')
ylabel('Mean NSE [dB]');
title({'Normalized Squared Error vs. $N_T, U$', ...
    ['$L = ', num2str(L) '; \delta_\Theta = ', num2str(delta_Theta, 3), '^{\circ}$']}, 'Interpreter', 'latex');


f = load('NSEs_NT_U_noisy_again.mat');
U_range = f.U_range;
N_T_range = f.N_T_range;

figure; hold on;
data = squeeze(f.NSEs_NT_NTheta_SS);
for U_ind = 1:length(U_range)
    plot(N_T_range, 10*log10(data(:, U_ind)), '-x');
end
legendString = ([repmat('$U = ', [3 1]), num2str([1 3 5].'), repmat(' $', [3 1])]);
legend(legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
xlabel('$N_T$', 'Interpreter','latex')
ylabel('Mean NSE [dB]');
title({'Normalized Squared Error vs. $N_T, U$', ...
    ['$L = ', num2str(L) '; \delta_\Theta = ', num2str(delta_Theta, 3), '^{\circ}$']}, 'Interpreter', 'latex');

% f = load("NSEs_NT_NTheta_noiseless.mat");
% N_Theta_range = f.N_Theta_range;
% N_T_range = f.N_T_range;
% 
% figure;
% hold on;
% data = squeeze(f.NSEs_NT_NTheta_SS);
% for N_Theta_ind = 1:length(N_Theta_range)
%     plot(N_T_range, 10*log10(data(N_Theta_ind, :)), '-x');
% end

