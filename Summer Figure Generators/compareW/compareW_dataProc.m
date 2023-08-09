f = load('compareW_1000_iters_NT_3.mat');

container = f.container;
SNR_dB_range = f.SNR_dB_range;

figure; hold on;
for i = 1:3
    curr_data = squeeze(container(i, :, :));
    mean_NSE_dB = 10*log10(mean(curr_data, 2));
    plot(SNR_dB_range, mean_NSE_dB, '-x');
end

xlabel('Per-receiver SNR $\gamma$','Interpreter', 'latex')
ylabel('Normalized Squared Error [dB]','Interpreter', 'latex')
legend({'Whole Data', 'Matched Filter', 'MMSE'}, 'Location','best', 'FontSize', 14);
grid on;

title({'Comparison of Imaging Post-Processing', '$N_T = 3; L =30, \delta_\Theta=.47^\circ$'}, 'interpreter', 'latex');