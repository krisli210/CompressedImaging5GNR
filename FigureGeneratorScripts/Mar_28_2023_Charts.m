close all
clear

% MNSE vs. N_T 
load("NSEs_N_T.mat");
figure;
plot(N_T_Range, mean(NSEs_N_T, 2), '-x') 
title('Normalized Squared Error vs. Time Slots')
xlabel('N_T');
ylabel('Mean NSE');
grid on;
xlim([1, inf]);

load("NSEs_SNRs.mat");
figure;
plot(SNR_range, mean(NSEs_SNRs, 2), '-x' , 'Color', 'r') 
title('Normalized Squared Error vs. Receive SNR')
xlabel('SNR [dB]');
ylabel('Mean NSE');
grid on;

load("NSEs_U_0SNR.mat");
figure;
plot(U_Range, mean(NSEs_U, 2), '-x', 'Color', 'm') 
title('Normalized Squared Error vs. U')
xlabel('U');
ylabel('Mean NSE');
grid on;