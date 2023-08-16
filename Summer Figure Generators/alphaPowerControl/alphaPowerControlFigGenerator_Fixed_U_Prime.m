close all
% clear

f = load('alphaPowerControl_vary_V_uniform_Scats_fix_U_prime_extra_Iters.mat');
f2 = load('alphaPowerControl_vary_V_uniform_Scats_fix_U_prime.mat')


alpha_range = f.alpha_range;
psnrs = f.psnrs;
sum_capacities = f.sum_capacities;
v_range = f.v_range;
u_range = f.u_range;

psnrs = zeros(length(v_range), length(alpha_range), 100+300);
sum_capacities = zeros(length(v_range), length(alpha_range), 100+300);

psnrs(:, :, 1:300) = f.psnrs;
psnrs(:, :, 301:400) = f2.psnrs; 

sum_capacities(:, :, 1:300) = f.sum_capacities;
sum_capacities(:, :, 301:400) = f2.sum_capacities;

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
legend({'$U=1$', '$U=4$', '$U=7$'}, 'Interpreter','latex','FontSize',14);
