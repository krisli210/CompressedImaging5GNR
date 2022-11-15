%% Comm-defined imaging of a MIMO scattering channel
close all
clear

rng(42);
% Begin with a PoC script utilizing randomized beamforming
% and 5GNR data packets formed with the 5GNR comm toolbox

% Scattering channel is captured by a MIMO monostatic RADAR
% Image is formed @ the RADAR via l1 optimization assuming sparse scenes
prm.CenterFreq = 28e9;

c = physconst('LightSpeed');
prm.PropagationSpeed = c;
prm.lam = c/prm.CenterFreq;

prm.BsPos = [0; 0; 0];
prm.BsArraySize = [8 8]; %BS Dimension
prm.NumBsElements = prod(prm.BsArraySize);
prm.BsAZlim = [-60 60];
prm.BsELlim = [-90 0];

prm.RxPos = [0; 0; 0];
prm.RxArraySize = [8 8];
prm.NumRxElements = prod(prm.RxArraySize);
prm.RxAZlim = [-60 60];
prm.RxELlim = [-90 0];

prm.NumUsers = 4;
prm.NumPackets = 100;
prm.Ns = 1; %number of symbols per packet
prm.M = 2; %modulation order
prm.K = 20^2; %give as square img 

%Arrays as uniform rectangular given in PA toolbox
BsArray = phased.URA(prm.BsArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));

RxArray = phased.URA(prm.RxArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement);

%Scatterer generation
nScat = 7;
rMin = 20; rMax = 80;
thetaMin = -60; %in Azimuth
thetaMax = 60;
[ScatPosPol, ScatPosCart] = generateScatPositions(nScat, rMin, rMax, thetaMin, thetaMax);
ScatCoeffs = complex(ones(1, nScat), ones(1, nScat)) ./ sqrt(2); %Unit reflector
%End Scatterer generation

% % % RADAR Channel Block
channelRADAR = phased.ScatteringMIMOChannel;
channelRADAR.PropagationSpeed = c;
channelRADAR.CarrierFrequency = prm.CenterFreq;
channelRADAR.SampleRate = 400e6;
channelRADAR.SimulateDirectPath = false;
channelRADAR.ChannelResponseOutputPort = true; %Gives delay as output
channelRADAR.Polarization = 'None'; %Dual polarization requires defining scattering functions in both horiz and vert directions

channelRADAR.TransmitArray = BsArray;
channelRADAR.TransmitArrayPosition = prm.BsPos; 
channelRADAR.ReceiveArray = RxArray;
channelRADAR.ReceiveArrayPosition = prm.RxPos; 
channelRADAR.ScattererSpecificationSource = 'Property';
channelRADAR.ScattererPosition = ScatPosCart;
channelRADAR.ScattererCoefficient = ScatCoeffs;
channelRADAR.MaximumDelaySource = 'Auto';

[~,trueH,tau] = channelRADAR(complex(randn(1,prm.NumBsElements), ...
    randn(1,prm.NumBsElements)));
trueH = sum(trueH, 3);
maxChDelay = ceil(max(tau)*channelRADAR.SampleRate);
% % % RADAR Channel Block

% [waveform, info] = genMU5GNRWaveform();
% txWaveform = repmat(waveform, [1 64]);

%tx signal construction
s = zeros(prm.NumUsers, prm.Ns * prm.NumPackets);
x = zeros(prm.NumBsElements, prm.Ns * prm.NumPackets);
qpskmod = comm.QPSKModulator('BitInput',true);   
for p = 1:prm.NumPackets
    F = complex(randn(prm.NumBsElements, prm.NumUsers), randn(prm.NumBsElements, prm.NumUsers));
    for u = 1:prm.NumUsers
        bits = randi(2, prm.M*prm.Ns, 1) - 1; 
        s(u, (p-1)*prm.Ns+1:p*prm.Ns) = qpskmod(bits);
    end
    x(:, (p-1)*prm.Ns+1:p*prm.Ns) = F * s(:, (p-1)*prm.Ns+1:p*prm.Ns);
end

% F = complex(randn(prm.NumBsElements, prm.NumUsers), randn(prm.NumBsElements, prm.NumUsers));
% x = [F*s zeros(prm.NumBsElements, maxChDelay)];
%end tx signal construction

[RefImg, Gamma, H_TX, H_RX, physH] = genRandomRefImage(prm, ScatPosPol, ScatPosCart, ScatCoeffs, rMax, thetaMin, thetaMax); 

% y = (channelRADAR(x.')).';
y = physH * x;
NNs = prm.NumRxElements * size(y, 2);
y_vec = reshape(y, [NNs 1]);

xKron = kron(x, eye(prm.NumRxElements));
H_combined = kr(H_TX, H_RX);
A = xKron.' * H_combined;

% x0 = .0001*ones(prm.K, 1);
x0 = complex(randn(prm.K, 1), randn(prm.K, 1));
% x0 = Gamma;
epsilon = 1e-4;
Gamma_hat = l1qc_logbarrier(x0, A, [], y_vec, epsilon);
% Gamma_hat = l1eq_pd(x0, A, 0, y_vec);
% Gamma_hat = l1dantzig_pd(x0, A, [], y_vec, epsilon)

figure;
subplot(1, 2, 1); imagesc(abs(RefImg).^2);
title('Reference Image'); xlabel('\theta'); ylabel('Range [m]');

img_hat = reshape(Gamma_hat, [sqrt(prm.K) sqrt(prm.K)]);
subplot(1, 2, 2); imagesc(abs(img_hat).^2);
title('Reconstructed Image'); xlabel('\theta'); ylabel('Range [m]');