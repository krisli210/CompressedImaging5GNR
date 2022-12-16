%% Comm-defined imaging of a MIMO scattering channel
close all
clear 

rng(423);
% Begin with a PoC script utilizing randomized beamforming
% and 5GNR data packets formed with the 5GNR comm toolbox

% Scattering channel is captured by a MIMO monostatic RADAR
% Image is formed @ the RADAR via l1 optimization assuming sparse scenes
prm.CenterFreq = 28e9;

c = physconst('LightSpeed');
prm.PropagationSpeed = c;
prm.lam = c/prm.CenterFreq;

prm.BsPos = [0; 0; 0];
prm.BsArraySize = [4 4]; %BS Dimension
prm.NumBsElements = prod(prm.BsArraySize);
prm.BsAZlim = [-60 60];
prm.BsELlim = [-90 0];

prm.RxPos = [0; 0; 0];
prm.RxArraySize = [4 4];
prm.NumRxElements = prod(prm.RxArraySize);
prm.RxAZlim = prm.BsAZlim;
prm.RxELlim = [-90 0];

prm.NumUsers = 4;
prm.NumPackets = 20;
prm.Ns = 1; %number of symbols per packet
prm.M = 2; %modulation order
prm.K = 16; %number of grid spots, i.e. dimension of the quantized azimuth profile

angles = prm.BsAZlim(1):(prm.BsAZlim(2)-prm.BsAZlim(1))/(prm.NumPackets-1):prm.BsAZlim(2);
%Arrays as uniform rectangular given in PA toolbox
% BsArray = phased.URA(prm.BsArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));
BsArray = phased.ULA(prm.NumBsElements, .5*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));

% RxArray = phased.URA(prm.RxArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement);
RxArray = phased.ULA(prm.NumRxElements, .5*prm.lam, 'Element', phased.IsotropicAntennaElement);

%Scatterer generation
nScat = 2;
rMin = 20; rMax = 80;
thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth
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
% trueH = sum(trueH, 3);
maxChDelay = ceil(max(tau)*channelRADAR.SampleRate);
% % % RADAR Channel Block

%DFT matrix definitions
A_tx = (1/sqrt(prm.NumBsElements)) .* exp(-1j * 2*pi * (0:prm.NumBsElements-1).' * (0:prm.NumBsElements-1) ./ prm.NumBsElements);
A_rx = (1/sqrt(prm.NumRxElements)) .* exp(-1j * 2*pi * (0:prm.NumRxElements-1).' * (0:prm.NumBsElements-1) ./ prm.NumRxElements);

%tx signal construction
s = zeros(prm.NumUsers, prm.Ns * prm.NumPackets);
x = zeros(prm.NumBsElements, prm.Ns * prm.NumPackets);
qpskmod = comm.QPSKModulator('BitInput',true);   
for p = 1:prm.NumPackets
%     F = fft(ones(prm.NumBsElements, prm.NumUsers));
%     F = complex(randn(prm.NumBsElements, prm.NumUsers), randn(prm.NumBsElements, prm.NumUsers));
%     F = steervec(getElementPosition(BsArray), angles(p));

    cols = randperm(prm.NumBsElements, prm.NumUsers); % DFT based precoding
    F = A_tx(:, cols);

    for u = 1:prm.NumUsers
        bits = randi(2, prm.M*prm.Ns, 1) - 1; 
        s(u, (p-1)*prm.Ns+1:p*prm.Ns) = qpskmod(bits);
    end
    x(:, (p-1)*prm.Ns+1:p*prm.Ns) = F * s(:, (p-1)*prm.Ns+1:p*prm.Ns);
    % x is NumBsElements x Ns
end
%end tx signal construction

% [RefImg, Gamma, H_TX, H_RX, physH] = genRandomRefImage(prm, ScatPosCart, ScatCoeffs, ...
%                                      getElementPosition(BsArray), getElementPosition(RxArray), ... 
%                                      rMax, thetaMin, thetaMax); 
 
BsSteer = phased.SteeringVector('SensorArray', BsArray);
RxSteer = phased.SteeringVector('SensorArray', RxArray);
[azProfile, H_TX, H_RX, physH] = genRandomAzProfile(prm, ...
                                                    getElementPosition(BsArray), ...
                                                    getElementPosition(RxArray), ...
                                                    thetaMin, thetaMax, ...
                                                    BsSteer, RxSteer);

% y = channelRADAR(x.').'; 
W = eye(size(H_RX, 1));
y = W * physH * x;
NNs = size(H_RX, 1) * size(y, 2);
y_vec = reshape(y, [NNs 1]);

Phi = kron(x.', W); % this is full rank necessarily
Psi = kr(H_TX, H_RX); % this is problematic
A = Phi * Psi;

% x0 = complex(randn(prm.K, 1), randn(prm.K, 1));
% epsilon = 1e-6;
% Gamma_hat = l1qc_logbarrier(x0, A, [], y_vec, epsilon);


% % % Native Solvers 

% sensingDict = sensingDictionary('CustomDictionary', A);
% [Gamma_hat, YI, I, R] = matchingPursuit(sensingDict, y_vec, maxIterations=100, Algorithm="OMP", maxerr={"L1", 1e-1});

% % % 


% figure;
% subplot(1, 2, 1); imagesc(abs(RefImg).^2);
% title('Reference Image'); xlabel('\theta'); ylabel('Range [m]');
% 
% img_hat = reshape(Gamma_hat, [sqrt(prm.K) sqrt(prm.K)]);
% subplot(1, 2, 2); imagesc(abs(img_hat).^2);
% title('Reconstructed Image'); xlabel('\theta'); ylabel('Range [m]');