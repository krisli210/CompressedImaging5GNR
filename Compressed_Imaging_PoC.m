%% Comm-defined imaging of a MIMO scattering channel
close all

rng(42);
% Begin with a PoC script utilizing randomized beamforming
% and 5GNR data packets formed with the 5GNR comm toolbox

% Scattering channel is captured by a MIMO monostatic RADAR
% Image is formed @ the RADAR via l1 optimization assuming sparse scenes
prm.CenterFreq = 28e9;

c = physconst('LightSpeed');
prm.PropagationSpeed = c;
lam = c/prm.CenterFreq;
prm.Lam = lam;

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
prm.Ns = 1000; %number of symbols
prm.M = 2; %modulation order
prm.K = 12^2; %give as square img 

%Arrays as uniform rectangular given in PA toolbox
BsArray = phased.URA(prm.BsArraySize, .5*lam, 'Element', phased.IsotropicAntennaElement('BackBaffled', true));

RxArray = phased.URA(prm.RxArraySize, .5*lam, 'Element', phased.IsotropicAntennaElement);

%Scatterer generation
nScat = 5;
rMin = 20; rMax = 80;
thetaMin = -60; %in Azimuth
thetaMax = 60;
[ScatPosPol, ScatPosCart] = generateScatPositions(nScat, rMin, rMax, thetaMin, thetaMax);
ScatCoeff = complex(ones(1, nScat), ones(1, nScat)) ./ sqrt(2); %Unit reflector
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
channelRADAR.ScattererCoefficient = ScatCoeff;
channelRADAR.MaximumDelaySource = 'Auto';

[~,~,tau] = channelRADAR(complex(randn(1,prm.NumBsElements), ...
    randn(1,prm.NumBsElements)));
maxChDelay = ceil(max(tau)*channelRADAR.SampleRate);
% % % RADAR Channel Block

% [waveform, info] = genMU5GNRWaveform();
% txWaveform = repmat(waveform, [1 64]);

%tx signal construction
s = zeros(prm.NumUsers, prm.Ns);
qpskmod = comm.QPSKModulator('BitInput',true);   
for i = 1:prm.NumUsers
    bits = randi(2, prm.M*prm.Ns, 1) - 1; 
    s(i, :) = qpskmod(bits);
end

[RefImg, H_TX, H_RX] = genRandomRefImage(prm.K, ScatPosPol, rMin, rMax, thetaMin, thetaMax); 
Gamma = reshape(img, [prm.K, 1]);

F = complex(randn(prm.NumBsElements, prm.NumUsers), randn(prm.NumBsElements, prm.NumUsers));
x = [F*s zeros(prm.NumBsElements, maxChDelay)];

y = (channelRADAR(x.')).';
NNs = prm.NumRxElements * size(y, 2);
y_vec = reshape(y, [NNs 1]);

xKron = kron(x.', eye(prm.NumRxElements));


function [ScatPosPol, ScatPosCart] = generateScatPositions(nScat, rMin, rMax, thetaMin, thetaMax)
    r = (rMax-rMin).*rand(1, nScat) + rMin;
    theta = (thetaMax-thetaMin).*rand(1, nScat) + thetaMin;
    
    ScatPosPol = [r; theta; zeros(1, nScat)];
    ScatPosCart = [r.*cosd(theta) ; r.*sind(theta) ; zeros(1, nScat)];
end

function gammas = solve_inverse(A, y)
    
end

