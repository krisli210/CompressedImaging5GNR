close all
% clear

% f = load('alphaPowerControl_vary_V_uniform_Scats_fix_U_prime_extra_Iters.mat');
% f2 = load('alphaPowerControl_vary_V_uniform_Scats_fix_U_prime.mat')
% 
% 
% alpha_range = f.alpha_range;
% psnrs = f.psnrs;
% sum_capacities = f.sum_capacities;
% v_range = f.v_range;
% u_range = f.u_range;
% 
% psnrs = zeros(length(v_range), length(alpha_range), 100+300);
% sum_capacities = zeros(length(v_range), length(alpha_range), 100+300);
% 
% psnrs(:, :, 1:300) = f.psnrs;
% psnrs(:, :, 301:400) = f2.psnrs; 
% 
% sum_capacities(:, :, 1:300) = f.sum_capacities;
% sum_capacities(:, :, 301:400) = f2.sum_capacities;
% 
% figure; hold on;
% for v_ind = 1:length(v_range)
%     psnrs_curr = squeeze(psnrs(v_ind, :, :));
%     sum_capacities_curr = squeeze(sum_capacities(v_ind, :, :));
%     plot(mean(psnrs_curr(1:19, :), 2), mean(sum_capacities_curr(1:19, :), 2), '-x');
% end
% 
% xlabel('Image PSNR [dB]','Interpreter','latex', 'FontSize',14);
% ylabel('Sum Spectral Efficiency [bps/Hz]', 'Interpreter','latex','FontSize',14);
% title('Virtual User Tradeoff Curve')
% grid on;
% legend({'$U=1$', '$U=4$', '$U=7$'}, 'Interpreter','latex','FontSize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% % clear
% 
f = load('alphaPowerControl_vary_V_uniform_Scats_fix_U_prime_rice_scene_v2.mat');


alpha_range = f.alpha_range;
psnrs = f.psnrs;
sum_capacities = f.sum_capacities;
v_range = f.v_range;
u_range = f.u_range;

% psnrs = zeros(length(v_range), length(alpha_range), 100+300);
% sum_capacities = zeros(length(v_range), length(alpha_range), 100+300);
% 
% psnrs(:, :, 1:300) = f.psnrs;
% psnrs(:, :, 301:400) = f2.psnrs; 

% sum_capacities(:, :, 1:300) = f.sum_capacities;
% sum_capacities(:, :, 301:400) = f2.sum_capacities;

figure; hold on;
for v_ind = 1:length(v_range)
    if (v_ind == 1)
        lin_str = '-x';
    elseif (v_ind == 2)
        lin_str = '-o';
    else
        lin_str = '-*';
    end
    psnrs_curr = squeeze(10.^(psnrs(v_ind, :, :)./10));
    sum_capacities_curr = squeeze(sum_capacities(v_ind, :, :));
    plot(10*log10(mean(psnrs_curr(1:19, :), 2)), mean(sum_capacities_curr(1:19, :), 2), lin_str);

    if (v_ind == 3)
        text(10*log10(mean(psnrs_curr(1, :), 2)), mean(sum_capacities_curr(1, :), 2), '$\alpha=1$', 'VerticalAlignment','bottom', 'HorizontalAlignment','left','Interpreter','latex','FontSize',14)
        text(10*log10(mean(psnrs_curr(19, :), 2)), mean(sum_capacities_curr(19, :), 2), '$\alpha=.1$', 'VerticalAlignment','bottom', 'HorizontalAlignment','right','Interpreter','latex','FontSize',14)
    end
end

xlabel('Image PSNR [dB]','Interpreter','latex', 'FontSize',14);
ylabel('$R_U$ [bps/Hz]', 'Interpreter','latex','FontSize',14);
ylim([-inf, 10]);
title('Tradeoff Curve for RICE Scene', 'FontSize', 14)
% title('Tradeoff Curve for DCM Scene', 'FontSize', 14)

grid on;
legend({'$U=1; V=7$', '$U=4; V=4$', '$U=7; V=1$'}, 'Interpreter','latex','FontSize',14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a verification of the power norms (renorm by sqrt(M)) and
% refactoring of the Montecarlo sim to meet Divyanshu recommendation to fix
% scenes rather than double averaging

% close all
% clear

f = load('alphaPowerControl_vary_V_uniform_Scats_fix_U_prime_uniform_scatterers_v4.mat');


alpha_range = f.alpha_range;
psnrs = f.psnrs;
sum_capacities = f.sum_capacities;
v_range = f.v_range;
u_range = f.u_range;

% psnrs = zeros(length(v_range), length(alpha_range), 100+300);
% sum_capacities = zeros(length(v_range), length(alpha_range), 100+300);
% 
% psnrs(:, :, 1:300) = f.psnrs;
% psnrs(:, :, 301:400) = f2.psnrs; 

% sum_capacities(:, :, 1:300) = f.sum_capacities;
% sum_capacities(:, :, 301:400) = f2.sum_capacities;

figure; hold on;
for v_ind = 1:length(v_range)
    if (v_ind == 1)
        lin_str = '-x';
    elseif (v_ind == 2)
        lin_str = '-o';
    else
        lin_str = '-*';
    end
    psnrs_curr = squeeze(psnrs(v_ind, :, :));
    sum_capacities_curr = squeeze(sum_capacities(v_ind, :, :));

    plot(mean(psnrs_curr(1:19, :), 2), mean(sum_capacities_curr(1:19, :), 2), lin_str);

    if (v_ind == 3)
        text(mean(psnrs_curr(1, :), 2), mean(sum_capacities_curr(1, :), 2), '$\alpha=1$', 'VerticalAlignment','bottom', 'HorizontalAlignment','left','Interpreter','latex','FontSize',14)
        text(mean(psnrs_curr(19, :), 2), mean(sum_capacities_curr(19, :), 2), '$\alpha=.1$', 'VerticalAlignment','bottom', 'HorizontalAlignment','left','Interpreter','latex','FontSize',14)
    end
end


xlabel('Image PSNR [dB]','Interpreter','latex', 'FontSize',14);
ylabel('$R_U$ [bps/Hz]', 'Interpreter','latex','FontSize',14);
ylim([-inf, 10]);
title('Tradeoff Curve for Random Scatterers','FontSize', 14)
grid on;
legend({'$U=1; V=7$', '$U=4; V=4$', '$U=7; V=1$'}, 'Interpreter','latex','FontSize',14);

