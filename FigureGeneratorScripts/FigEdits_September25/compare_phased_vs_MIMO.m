close all
f = load("compare_phased_vs_MIMO_U_SNR_N_T_1.mat");

data = squeeze(f.wholeData);
delta_Theta = 120/f.prm.N_theta;

U_range = f.U_range;
SNR_range = f.SNR_dB_range;

figure; hold on;
for U_ind = 1:length(U_range)
    if(U_ind==1)
        lin_str = 'r-x';
    elseif(U_ind==2)
        lin_str = 'g-o';
    else
        lin_str = 'b-*';
    end
    plot(SNR_range(1:9), 10*log10(mean(data(U_ind, 1:9, :), 3)), lin_str);
end
grid on;

legendString = ([repmat('$U = ', [3 1]), num2str(U_range.'), repmat('$', [3 1])]);
legend(legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 24);
xlabel('Receive SNR [dB]', 'Interpreter','latex', 'FontSize',24)
ylabel('Mean NSE [dB]', 'Interpreter','latex', 'FontSize',24);
% title({'Image NSE vs. $SNR, U$', ...
%     ['$N_T = ', num2str(1), ';L = ', num2str(50) '; \delta_\Theta = ', num2str(delta_Theta, 3), '^{\circ}$']}, 'Interpreter', 'latex');
