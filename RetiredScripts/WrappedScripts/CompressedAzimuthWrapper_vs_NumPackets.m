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

KRange = [64 128 256 512];


nPerNumPackets = 1000;
packetRange = 1:20;

for KIndex = 1:length(KRange)
    prm.K = KRange(KIndex);

    Pd = zeros(1, length(packetRange));
    Pf = zeros(1, length(packetRange));
    for packetIndex = 1:length(packetRange)
        prm.NumPackets = packetRange(packetIndex);
        prm.Ns = ceil(40/prm.NumPackets);
        PdPerPacket = 0;
        PfPerPacket = 0;
        for n = 1:nPerNumPackets
            [detect, falseAlarm] = CompressedAzimuthWrapped(prm);
            detectProbability = sum(detect) / prm.NumTargets;
            
            PdPerPacket = PdPerPacket + detectProbability;
            PfPerPacket = PfPerPacket + falseAlarm;
        end
        Pd(packetIndex) = PdPerPacket / nPerNumPackets;
        Pf(packetIndex) = PfPerPacket / nPerNumPackets;
        
    end
    save(['data/azProfile/Pd_', 'K_', num2str(KRange(KIndex)), 'NPacket_', ...
            num2str(packetRange(packetIndex))], 'Pd')

    save(['data/azProfile/Pf_', 'K_', num2str(KRange(KIndex)), 'NPacket_', ...
            num2str(packetRange(packetIndex))], 'Pf')
end
