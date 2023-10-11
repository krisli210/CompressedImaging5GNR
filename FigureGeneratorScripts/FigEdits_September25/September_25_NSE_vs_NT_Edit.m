f = load('NSEs_N_T_U_noiseless.mat');
U_range = f.U_range;
N_T_range = f.N_T_range;

figure; hold on;
data = squeeze(f.NSEs_L_U_SS);
L = f.prm.L;
delta_Theta  = 120 / f.N_Theta_range;
for U_ind = 1:length(U_range)

    if(U_ind==1)
        lin_str = 'r-x';
    elseif(U_ind==2)
        lin_str = 'g-o';
    else
        lin_str = 'b-*';
    end
    plot(N_T_range, 10*log10(data(:, U_ind)), lin_str);
end
legendString = ([repmat('$U = ', [3 1]), num2str([1 3 5].'), repmat(' $', [3 1])]);
legend(legendString, 'location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
xlabel('Time Slots', 'Interpreter','latex', 'FontSize',14)
ylabel('Mean NSE [dB]', 'Interpreter', 'latex', 'FontSize', 14);
title({'Image NSE vs. $N_T$', ...
    ['$\delta_\Theta = ', num2str(delta_Theta, 3), '^{\circ}$']}, ...
    'FontSize' ,14, 'Interpreter', 'latex');