close all
clear
% 
f_interference = load('C:\Users\krisl\Desktop\Summer2023\CompressedImaging5GNR\Summer Figure Generators\alphaPowerControl\alphaPowerControl_vary_V_RICEscene_Interference_v2.mat');
% f_interference = load('C:\Users\krisl\Desktop\Summer2023\CompressedImaging5GNR\Summer Figure Generators\alphaPowerControl\alphaPowerControl_vary_V_UniformScene_Interference_v2.mat');

% f_noInterference = load("C:\Users\krisl\Desktop\Summer2023\CompressedImaging5GNR\Summer Figure Generators\alphaPowerControl\alphaPowerControl_vary_V_RICEscene_noInterference_v2.mat");

alpha_range = f_interference.alpha_range;
v_range = f_interference.v_range;
u_range = f_interference.u_range;

psnrs_interference = f_interference.psnrs;
sum_capacities = f_interference.sum_capacities;

% psnrs_noInterference = f_noInterference.psnrs;
% sum_capacities_noInterference = f_noInterference.sum_capacities;

figure; hold on;
for v_ind = 1:length(v_range)
    if(v_ind==1)
        lin_str = 'r-x';
    elseif (v_ind == 2)
        lin_str = 'g-o';
    else
        lin_str = 'b-*';
    end

    psnrs_curr_interference = squeeze(10.^(psnrs_interference(v_ind, :, :)./10));
    sum_capacities_curr_interference = squeeze(sum_capacities(v_ind, :, :));
    plot(10*log10(mean(psnrs_curr_interference(2:19, :), 2)), mean(sum_capacities_curr_interference(2:19, :), 2), ...
        lin_str);

    if (v_ind == 2)
        text(10*log10(mean(psnrs_curr_interference(2, :), 2)), mean(sum_capacities_curr_interference(2, :), 2), '$\alpha=.95$', 'VerticalAlignment','bottom', 'HorizontalAlignment','left','Interpreter','latex','FontSize',14)
    end
    if (v_ind == 2)
        text(10*log10(mean(psnrs_curr_interference(19, :), 2)), mean(sum_capacities_curr_interference(19, :), 2), '$\alpha=.1$', 'VerticalAlignment','bottom', 'HorizontalAlignment','right','Interpreter','latex','FontSize',14)
    end
end

% for v_ind = 1:length(v_range)
%     if(v_ind==1)
%         lin_str = 'r--x';
%     elseif (v_ind == 2)
%         lin_str = 'g--o';
%     else
%         lin_str = 'b--*';
%     end
%     psnrs_curr_noInterference = squeeze(10.^(psnrs_noInterference(v_ind, :, :)./10));
%     sum_capacities_curr_noInterference = squeeze(sum_capacities_noInterference(v_ind, :, :));
%     plot(10*log10(mean(psnrs_curr_noInterference(1:19, :), 2)), mean(sum_capacities_curr_noInterference(1:19, :), 2), ...
%         lin_str);
% end 

xlabel('Image PSNR [dB]','Interpreter','latex', 'FontSize',14);
ylabel('$R_U$ [bps/Hz]', 'Interpreter','latex','FontSize',14);
% ylim([-inf, 10]);
title('Tradeoff Curve for RICE Scene', 'FontSize', 14)
% title('Tradeoff Curve for DCM Scene', 'FontSize', 14)

grid on;
legend({'$U=1; V=7$', '$U=4; V=4$', '$U=7; V=1$'}, 'Interpreter','latex','FontSize',14, 'Location','best');