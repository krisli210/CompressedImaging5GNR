close all
clear
% 
% f_interference = load('C:\Users\krisl\Desktop\Summer2023\CompressedImaging5GNR\Summer Figure Generators\alphaPowerControl\alphaPowerControl_vary_V_RICEscene_Interference_v2.mat');
% f_interference = load('alphaPowerControl_vary_V_RICEscene_Interference_v3_0dBGain_fullData_n100.mat');
f_interference = load('alphaPowerControl_vary_V_RICEscene_Interference_v3_0dBGain_v4.mat');


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
x = [];
y = [];
for v_ind = [2, 3, 4]
    if(v_ind==2)
        lin_str = 'r--x';
    elseif (v_ind == 3)
        lin_str = 'g--o';
    else
        lin_str = 'b--*';
    end

    psnrs_curr_interference = squeeze(10.^(psnrs_interference(v_ind, 2:19, :)./10));
    sum_capacities_curr_interference = squeeze(sum_capacities(v_ind, 2:19, :));
    % plot(10*log10(mean(psnrs_curr_interference(:, :), 2)), mean(sum_capacities_curr_interference(:, :), 2), ...
    %     lin_str);
    plot(10*log10(mean(psnrs_curr_interference(:, :), 2)), mean(sum_capacities_curr_interference(:, :), 2), lin_str);
end

for v_ind = 2:length(v_range)-1
    % if(v_ind==1)
    %     lin_str = 'r-x';
    % elseif (v_ind == 2)
    %     lin_str = 'g-o';
    % else
    %     lin_str = 'b-*';
    % end

    psnrs_curr_interference = squeeze(10.^(psnrs_interference(v_ind, 2:19, :)./10));
    sum_capacities_curr_interference = squeeze(sum_capacities(v_ind, 2:19, :));
    
    x = [x; 10*log10(mean(psnrs_curr_interference, 2))];
    y = [y; mean(sum_capacities_curr_interference, 2)];
end
x = [x; min(x, [], 'all')];
y = [y; min(y, [], 'all')];
convexHull = convhull(x, y);
plot(x(convexHull), y(convexHull), 'Color', [0 0 0]);

xlabel('Image PSNR [dB]','Interpreter','latex', 'FontSize',14);
ylabel('$R_U$ [bps/Hz]', 'Interpreter','latex','FontSize',14);
xlim([min(x, [], 'all'), inf]);
ylim([min(y, [], 'all'), inf]);
% ylim([-inf, 10]);
% title({'Tradeoff Curve for RICE Scene', '$U+V = 8$'}, 'FontSize', 14, 'Interpreter', 'latex')

grid on;
legend({'$U=1; V=7$', '$U=2; V=6$', '$U=3; V=5$', 'Convex Hull'}, 'Interpreter','latex','FontSize',14, 'Location','best');