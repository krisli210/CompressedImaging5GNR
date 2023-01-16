function [detect, falseAlarm] = CompressedAzimuthWrapped(prm)
% Image is formed @ the RADAR via l1 optimization assuming sparse scenes
%Arrays as uniform rectangular given in PA toolbox
thetaMin = prm.BsAZlim(1); thetaMax = prm.BsAZlim(2); %in Azimuth

[azProfile, H_TX, H_RX, physH] = genRandomAzProfile(prm, ...
                                                    thetaMin, thetaMax ...
                                                    );

%tx signal construction
x = constructTxSignal(prm, H_TX);
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

[detect, falseAlarm] = measureDetection(azProfile, z_hat);

% % % 



end