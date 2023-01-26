close all

Kstr = string([repmat('K = ', 3, 1), num2str([64; 128; 256])]);
KRange =[64 128 256];

userRange = 1:16;

figure; hold on; grid on;
for KIndex = 1:length(KRange)

    Pd = load(['data/azProfile/Pd_', 'K_', num2str(KRange(KIndex)), 'U', ...
        num2str(16)], 'Pd').Pd;

    plot(userRange, Pd, '-x');
end
title('Probability of Detection vs. Users Per Channel')
xlabel('U');
ylabel('P_D')
legend(Kstr, 'Location', 'Southeast');

figure; hold on; grid on;
for KIndex = 1:length(KRange)
    
    Pf = load(['data/azProfile/Pf_', 'K_', num2str(KRange(KIndex)), 'U', ...
        num2str(16)], 'Pf').Pf;
    plot(userRange, Pf, '-x');
end
xlabel('U')
ylabel('P_F')
title('Probability of False Alarm vs. Users Per Channel')
legend(Kstr, 'Location', 'Northeast');