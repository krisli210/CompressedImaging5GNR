close all

Kstr = string([repmat('K = ', 4, 1), num2str([64; 128; 256; 512])]);
KRange =[64 128 256 512];

packetRange = 1:20;

figure; hold on; grid on;
for KIndex = 1:length(KRange)

    Pd = load(['data/azProfile/FixedNNs_vs_NPacket/Pd_', 'K_', num2str(KRange(KIndex)), 'NPacket_', ...
        num2str(20)], 'Pd').Pd;

    plot(packetRange, Pd, '-x');
end
title('Probability of Detection vs. Distinct Spatial Samples')
xlabel('Number of Precoders Measured');
ylabel('P_D')
legend(Kstr, 'Location', 'Southeast');

figure; hold on; grid on;
for KIndex = 1:length(KRange)
    
    Pf = load(['data/azProfile/FixedNNs_vs_NPacket/Pf_', 'K_', num2str(KRange(KIndex)), 'NPacket_', ...
        num2str(20)], 'Pf').Pf;
    plot(packetRange, Pf, '-x');
end
title('Probability of False Alarm vs. Distinct Spatial Samples')
xlabel('Number of Precoders Measured')
ylabel('P_F')
legend(Kstr, 'Location', 'Northeast');