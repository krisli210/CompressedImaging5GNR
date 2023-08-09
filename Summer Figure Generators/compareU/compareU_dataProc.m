f = load('compareU_1000_iters_NT_1.mat');

container = f.container;
SNR_dB_range = f.SNR_dB_range;
U_range = f.U_range;


figure; hold on;
for u = 1:3
    curr_data = squeeze(container(u, :, :));
    mean_NSE_dB = 10*log10(mean(curr_data, 2));
    plot(SNR_dB_range, mean_NSE_dB, '-x');
end

xlabel('Per-receiver SNR $\gamma$','Interpreter', 'latex')
ylabel('Normalized Squared Error [dB]','Interpreter', 'latex')
legend({'$U = 1$', '$U=6$', '$U=11$'}, 'Location','best', 'FontSize', 14, 'Interpreter','latex');
grid on;

title({'Comparison of Imaging Post-Processing', '$N_T = 1; L =30, \delta_\Theta=.47^\circ$'}, 'interpreter', 'latex');