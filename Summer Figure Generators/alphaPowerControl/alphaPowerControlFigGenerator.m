close all
% clear

f = load('alphaPowerControl_vary_V_uniform_Scats.mat');

alpha_range = f.alpha_range;
psnrs = f.psnrs;
sum_capacities = f.sum_capacities;
v_range = f.v_range;

figure; hold on;
for v_ind = 1:length(v_range)
    psnrs_curr = squeeze(psnrs(v_ind, :, :));
    sum_capacities_curr = squeeze(sum_capacities(v_ind, :, :));
    plot(mean(psnrs_curr(1:19, :), 2), mean(sum_capacities_curr(1:19, :), 2), '-x');
end

xlabel('Image PSNR [dB]','Interpreter','latex', 'FontSize',14);
ylabel('Sum Spectral Efficiency [bps/Hz]', 'Interpreter','latex','FontSize',14);
title('Virtual User Tradeoff Curve')
grid on;
legend({'$V=1$', '$V=3$', '$V=5$'}, 'Interpreter','latex','FontSize',14);
