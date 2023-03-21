clear 
close all

rng(42);

prm.CenterFreq = 28e9;

c = physconst('LightSpeed');
prm.PropagationSpeed = c;
prm.lam = c/prm.CenterFreq;

prm.BsPos = [0; 0; 0];
prm.BsArraySize = 16; %BS Dimension
prm.NumBsElements = prod(prm.BsArraySize);
prm.DeltaT = .5; % Element spacing normalized by wavelength
prm.BsAZlim = [-60 60];
prm.BsELlim = [-90 0];

prm.RxPos = [0; 0; 0];
prm.RxArraySize = 16;
prm.DeltaR = prm.NumBsElements * .5; % Set for virtual array 
prm.NumRxElements = prod(prm.RxArraySize);
prm.RxAZlim = prm.BsAZlim;
prm.RxELlim = [-90 0];

prm.NumUsers = 4;
prm.NumPackets = 3;
prm.Ns = 10; %number of symbols per packet
prm.M = 2; %modulation order
prm.K = 256; %number of grid spots, i.e. dimension of the quantized azimuth profile
prm.NumTargets = 3;

KRange = 64:16:384;
packetRange = 1:2:19;
userRange = 1:2:9;
NTargetRange = 1:2:9;

nPerNumPackets = 1000;


PdGrid = zeros([length(KRange), length(packetRange), length(userRange), length(NTargetRange)]);
PfGrid = zeros([length(KRange), length(packetRange), length(userRange), length(NTargetRange)]);

for KIndex = 1:length(KRange)
    prm.K = KRange(KIndex);

    for packetIndex = 1:length(packetRange)
        prm.NumPackets = packetRange(packetIndex);

        for userIndex = 1:length(userRange)
            prm.NumUsers = userRange(userIndex);

            for NTargetIndex = 1:length(NTargetRange)
            prm.NumTargets = NTargetRange(NTargetIndex);

            disp(['On: K - ', num2str(prm.K), ', NPacket - ', num2str(prm.NumPackets), ...
                ', NumUsers - ', num2str(prm.NumUsers), ', NumTargets - ', num2str(prm.NumTargets)]);

            PdPerWrap = 0;
            PfPerWrap = 0;
                for n = 1:nPerNumPackets
                    [detect, falseAlarm] = CompressedAzimuthWrapped(prm);
                    detectProbability = sum(detect) / prm.NumTargets;
                    
                    PdPerWrap = PdPerWrap + detectProbability;
                    PfPerWrap = PfPerWrap + falseAlarm;
                end
            PdGrid(KIndex, packetIndex, userIndex, NTargetIndex) = PdPerWrap / nPerNumPackets;
            PfGrid(KIndex, packetIndex, userIndex, NTargetIndex) = PfPerWrap / nPerNumPackets;
            end
        end
    end
end
save('data/azProfile/OMP_GridSearch/PdGrid', ...
            'PdGrid', 'KRange', 'packetRange', 'userRange', 'NTargetRange');

save('data/azProfile/OMP_GridSearch/PfGrid', ...
            'PfGrid', 'KRange', 'packetRange', 'userRange', 'NTargetRange');
