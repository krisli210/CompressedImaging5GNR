%% Comm-defined imaging of a MIMO scattering channel
close all
% clear 

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
prm.NumPackets = 4;
prm.Ns = 1; %number of symbols per packet
prm.M = 2; %modulation order
prm.K = 64; %number of grid spots, i.e. dimension of the quantized azimuth profile

angles = prm.BsAZlim(1):(prm.BsAZlim(2)-prm.BsAZlim(1))/(prm.NumPackets-1):prm.BsAZlim(2);
%Arrays as uniform rectangular given in PA toolbox
% BsArray = phased.URA(prm.BsArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));
BsArray = phased.ULA(prm.NumBsElements, prm.DeltaT*prm.lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));

% RxArray = phased.URA(prm.RxArraySize, .5*prm.lam, 'Element', phased.IsotropicAntennaElement);
RxArray = phased.ULA(prm.NumRxElements, prm.DeltaR*prm.lam, 'Element', phased.IsotropicAntennaElement);

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

% [RefImg, z, H_TX, H_RX, physH] = genRandomRefImage(prm, ScatPosCart, ScatCoeffs, ...
%                                      getElementPosition(BsArray), getElementPosition(RxArray), ... 
%                                      rMax, thetaMin, thetaMax); 

[azProfile, H_TX, H_RX, physH] = genRandomAzProfile(prm, ...
                                                    thetaMin, thetaMax ...
                                                    );

% figure;
% pattern(BsArray,prm.CenterFreq, Weights=H_TX, EL=0)

%tx signal construction
x = constructTxSignal(prm, H_TX);

% y = channelRADAR(x.').'; 
W = eye(size(H_RX, 1));
y = W * physH * x;
NNs = size(H_RX, 1) * size(y, 2);
y_vec = reshape(y, [NNs 1]);

Phi = kron(x.', W); % this is full rank necessarily
Psi = kr(H_TX, H_RX); % 
% Psi = kron(H_TX', H_RX);
A = Phi * Psi;

% % % Native Solvers 

sensingDict = sensingDictionary('CustomDictionary', A);
[z_hat, YI, I, R] = matchingPursuit(sensingDict, y_vec, maxIterations=10, Algorithm="OMP", maxerr={"L1", 1e-4});

A_sub = Phi * Psi(:, I); % Solve magnitude posthence via direct linsolve against estimated support
mags = linsolve(A_sub, y_vec);
z_hat(I) = mags;

% if  (find(z_hat) == find(azProfile))
%     disp('Correct AoA Estimation');
% end
figure; hold on; 
stem(-60:120/(prm.K-1):60, abs(azProfile));
stem(-60:120/(prm.K-1):60, abs(z_hat), '--x');

legend({'True', 'Estimate'})
% Magnitude issue must fall to the MP algorithm since everythign else is just thrown into A 

% % % 




% figure;
% subplot(1, 2, 1); imagesc(abs(RefImg).^2);
% title('Reference Image'); xlabel('\theta'); ylabel('Range [m]');
% 
% img_hat = reshape(z_hat, [sqrt(prm.K) sqrt(prm.K)]);
% subplot(1, 2, 2); imagesc(abs(img_hat).^2);
% title('Reconstructed Image'); xlabel('\theta'); ylabel('Range [m]');