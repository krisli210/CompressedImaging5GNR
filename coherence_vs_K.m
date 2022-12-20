close all
clear 

figure; hold on;

arr_dim = 16:16:64;
K = 8:256;
arr_ind = 1;
for arr = arr_dim
    coh_phi = zeros(length(K), 1);
    k_ind = 1;
    for k = K
        prm.CenterFreq = 28e9;
        
        c = physconst('LightSpeed');
        prm.PropagationSpeed = c;
        prm.lam = c/prm.CenterFreq;
        
        prm.BsPos = [0; 0; 0];
        prm.NumBsElements = arr;
        prm.BsAZlim = [-60 60];
        prm.BsELlim = [-90 0];
        
        prm.RxPos = [0; 0; 0];
        prm.NumRxElements = arr;
        prm.RxAZlim = prm.BsAZlim;
        prm.RxELlim = [-90 0];
        
        prm.NumUsers = 4;
        prm.NumPackets = 20;
        prm.Ns = 1; %number of symbols per packet
        prm.M = 2; %modulation order
        prm.K = k; %number of grid spots, i.e. dimension of the quantized azimuth profile
        
        thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
        
        angles = prm.BsAZlim(1):(prm.BsAZlim(2)-prm.BsAZlim(1))/(prm.NumPackets-1):prm.BsAZlim(2);
        %Arrays as uniform rectangular given in PA toolbox
        % BsArray = phased.URA(prm.BsArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));
        BsArray = phased.ULA(prm.NumBsElements, .5*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));
        
        % RxArray = phased.URA(prm.RxArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement);
        RxArray = phased.ULA(prm.NumRxElements, .5*prm.lam, 'Element', phased.IsotropicAntennaElement);
        
        
        BsSteer = phased.SteeringVector('SensorArray', BsArray);
        RxSteer = phased.SteeringVector('SensorArray', RxArray);
        [azProfile, H_TX, H_RX, physH] = genRandomAzProfile(prm, ...
                                                            getElementPosition(BsArray), ...
                                                            getElementPosition(RxArray), ...
                                                            thetaMin, thetaMax, ...
                                                            BsSteer, RxSteer);
        coh_phi(k_ind) = mutual_coherence(kr(H_TX, H_RX));
        k_ind = k_ind + 1;
    end
    plot(K, coh_phi)
    arr_ind = arr_ind + 1;
end
title('Mutual Coherence of Dictionary vs. K');
xlabel('K');
ylabel('\mu(\Psi)');
leg_string = string([repmat('M = N = ', 4, 1), num2str([16; 32; 48; 64])]);
